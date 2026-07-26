"""Post-hoc QA for GutIBM AWS Batch array outputs.

Downloads per-index ``output.h5.gz`` (and optional ``input.json``), gunzips,
reads summary scalars, and checks that distinct seeds produced distinct
outcomes. This is the Stage 3 gate template: prove the array actually delivered
usable, seed-sensitive artifacts before trusting campaign science.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .hdf5_reader import GutIBMData
from .path_utils import (
    PathValidationError,
    prepare_output_directory,
    prepare_output_file,
    validate_path_syntax,
)

AWS = "aws"
INPUT_FILE_NAME = "input.json"
OUTPUT_FILE_NAME = "output.h5.gz"
H5_FILE_NAME = "output.h5"

AwsRunner = Callable[[list[str]], str]
S3Download = Callable[[str, Path], None]


@dataclass(frozen=True)
class IndexSummary:
    index: int
    seed: int | None
    final_agents: int
    final_time_s: float
    boundary_deaths: int
    washout_deaths: int
    colicin_kills: int
    mean_carbon: float
    fingerprint: str
    has_agents_layer: bool
    steps: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class ArrayQaReport:
    rows: list[IndexSummary]
    distinct_fingerprints: bool
    distinct_seeds: bool
    warnings: list[str] = field(default_factory=list)


def _default_aws_run(argv: list[str]) -> str:
    result = subprocess.run(argv, capture_output=True, text=True, check=True)
    return result.stdout


def _default_s3_download(s3_uri: str, dest: Path) -> None:
    prepare_output_file(dest)
    result = subprocess.run(
        [AWS, "s3", "cp", s3_uri, str(dest)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise FileNotFoundError(f"{s3_uri}: {result.stderr.strip()}")


def _strip_s3_prefix(prefix: str) -> str:
    text = prefix.strip()
    if not text.startswith("s3://"):
        raise ValueError(f"expected s3:// URI, got {prefix!r}")
    return text.rstrip("/")


def list_array_indices(
    output_prefix: str,
    *,
    aws_run: AwsRunner = _default_aws_run,
) -> list[int]:
    """List numeric child indexes under an array output prefix."""
    prefix = _strip_s3_prefix(output_prefix) + "/"
    listing = aws_run([AWS, "s3", "ls", prefix])
    indices: list[int] = []
    for line in listing.splitlines():
        token = line.strip().rstrip("/")
        if not token:
            continue
        # ``aws s3 ls s3://bucket/prefix/`` yields ``PRE 0/`` lines.
        name = token.split()[-1]
        if name.isdigit():
            indices.append(int(name))
    return sorted(set(indices))


def _gunzip_h5(gz_path: Path, h5_path: Path) -> None:
    prepare_output_file(h5_path)
    with gzip.open(gz_path, "rb") as src, h5_path.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def _event_count(summary: dict[str, Any], name: str) -> int:
    events = summary.get("events") or {}
    value = events.get(name, 0)
    return int(value)


def _fingerprint(summary: dict[str, Any], seed: int | None) -> str:
    payload = {
        "seed": seed,
        "time": summary.get("time"),
        "step": summary.get("step"),
        "n_total": summary.get("n_total", summary.get("num_agents")),
        "events": summary.get("events", {}),
        "chem": summary.get("chem", {}),
    }
    blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:16]


def summarize_hdf5(path: Path, *, index: int, seed: int | None) -> IndexSummary:
    """Build an IndexSummary from a local HDF5 file."""
    validate_path_syntax(path)
    with GutIBMData(path) as data:
        steps = data.steps
        if not steps:
            raise ValueError(f"{path}: no summary steps")
        final = data.get_summary(steps[-1])
        has_agents = False
        if data._file is not None and "agents" in data._file:
            has_agents = any(
                name.startswith("step_") for name in data._file["agents"]
            )
    return IndexSummary(
        index=index,
        seed=seed,
        final_agents=int(final.get("n_total", final.get("num_agents", 0))),
        final_time_s=float(final.get("time", 0.0)),
        boundary_deaths=_event_count(final, "boundary_deaths"),
        washout_deaths=_event_count(final, "washout_deaths"),
        colicin_kills=_event_count(final, "colicin_kills"),
        mean_carbon=float((final.get("chem") or {}).get("mean_carbon", 0.0)),
        fingerprint=_fingerprint(final, seed),
        has_agents_layer=has_agents,
        steps=list(steps),
    )


def _read_seed(input_path: Path) -> int | None:
    cfg = json.loads(input_path.read_text(encoding="utf-8"))
    seed = cfg.get("seed")
    return int(seed) if seed is not None else None


def qa_array(
    *,
    output_prefix: str,
    input_prefix: str | None = None,
    indices: list[int] | None = None,
    work_dir: Path,
    require_distinct: bool = True,
    aws_run: AwsRunner = _default_aws_run,
    s3_download: S3Download = _default_s3_download,
) -> ArrayQaReport:
    """Download array outputs and compare summary fingerprints."""
    out_root = _strip_s3_prefix(output_prefix)
    in_root = _strip_s3_prefix(input_prefix) if input_prefix else None
    work = prepare_output_directory(work_dir)

    resolved = indices if indices is not None else list_array_indices(
        out_root, aws_run=aws_run
    )
    if len(resolved) < 1:
        raise ValueError(f"no array indexes found under {out_root}")

    rows: list[IndexSummary] = []
    warnings: list[str] = []
    for index in resolved:
        child = prepare_output_directory(work / str(index))
        gz_path = child / OUTPUT_FILE_NAME
        h5_path = child / H5_FILE_NAME
        s3_download(f"{out_root}/{index}/{OUTPUT_FILE_NAME}", gz_path)
        _gunzip_h5(gz_path, h5_path)

        seed: int | None = None
        if in_root is not None:
            input_path = child / INPUT_FILE_NAME
            s3_download(f"{in_root}/{index}/{INPUT_FILE_NAME}", input_path)
            seed = _read_seed(input_path)

        row = summarize_hdf5(h5_path, index=index, seed=seed)
        if not row.has_agents_layer:
            warnings.append(
                f"index {index}: summary-only HDF5 (no agents layer) — "
                "spatial analysis/validation tools cannot run yet"
            )
        rows.append(row)

    fingerprints = {row.fingerprint for row in rows}
    seeds = {row.seed for row in rows if row.seed is not None}
    distinct_fp = len(fingerprints) == len(rows)
    distinct_seeds = len(seeds) == len(rows) if seeds else False

    if require_distinct and len(rows) >= 2 and not distinct_fp:
        raise AssertionError(
            "array outputs are not distinct: summary fingerprints collide "
            f"({sorted(fingerprints)})"
        )
    if require_distinct and in_root is not None and len(rows) >= 2 and not distinct_seeds:
        warnings.append("input seeds are missing or not all distinct")

    return ArrayQaReport(
        rows=rows,
        distinct_fingerprints=distinct_fp,
        distinct_seeds=distinct_seeds,
        warnings=warnings,
    )


def format_report(report: ArrayQaReport) -> str:
    """Human-readable table for CLI / docs paste."""
    lines = [
        "index  seed   agents  t_final  boundary  washout  colicin  fingerprint        agents_h5",
        "-----  -----  ------  -------  --------  -------  -------  -----------------  ---------",
    ]
    for row in report.rows:
        seed = "-" if row.seed is None else str(row.seed)
        lines.append(
            f"{row.index:<5d}  {seed:<5s}  {row.final_agents:<6d}  "
            f"{row.final_time_s:<7.0f}  {row.boundary_deaths:<8d}  "
            f"{row.washout_deaths:<7d}  {row.colicin_kills:<7d}  "
            f"{row.fingerprint:<17s}  {'yes' if row.has_agents_layer else 'no'}"
        )
    lines.append("")
    lines.append(
        f"distinct_fingerprints={report.distinct_fingerprints}  "
        f"distinct_seeds={report.distinct_seeds}"
    )
    for warning in report.warnings:
        lines.append(f"warning: {warning}")
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "QA GutIBM Batch array outputs: download, gunzip, compare summary "
            "fingerprints across indexes (seed sensitivity gate)."
        )
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="S3 prefix with per-index outputs, e.g. s3://bucket/.../out",
    )
    parser.add_argument(
        "--input-prefix",
        default=None,
        help="Optional S3 prefix with per-index input.json (for seed join)",
    )
    parser.add_argument(
        "--indices",
        default=None,
        help="Comma-separated indexes (default: list children under output prefix)",
    )
    parser.add_argument(
        "--work-dir",
        default="batch_qa_work",
        help="Local directory for downloads (default: batch_qa_work)",
    )
    parser.add_argument(
        "--allow-identical",
        action="store_true",
        help="Do not fail when fingerprints collide (report only)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    indices = None
    if args.indices:
        indices = [int(part.strip()) for part in args.indices.split(",") if part.strip()]
    try:
        report = qa_array(
            output_prefix=args.output_prefix,
            input_prefix=args.input_prefix,
            indices=indices,
            work_dir=Path(args.work_dir),
            require_distinct=not args.allow_identical,
        )
    except (AssertionError, FileNotFoundError, PathValidationError, ValueError, OSError) as exc:
        print(f"aws batch qa error: {exc}", file=sys.stderr)
        return 2
    print(format_report(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
