"""CLI for iterative GutIBM checkpoint analysis.

Modes:
  --summarize   Stream checkpoints → summary.csv (+ JSON metadata)
  --snapshot    One checkpoint → spatial PNGs
  --timeseries  summary.csv → trend PNGs
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

from .checkpoint_plots import write_snapshot_figures, write_timeseries_figures
from .checkpoint_scan import discover_checkpoints, validate_checkpoint_dir
from .checkpoint_summary import extract_many
from .path_utils import (
    PathValidationError,
    prepare_output_directory,
    prepare_output_file,
    validate_input_path,
    validate_path_syntax,
)

MODE_SUMMARIZE = "summarize"
MODE_SNAPSHOT = "snapshot"
MODE_TIMESERIES = "timeseries"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gut-ibm-analyze",
        description="Analyze GutIBM Spec-4 HDF5 checkpoints one file at a time.",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--summarize",
        action="store_const",
        const=MODE_SUMMARIZE,
        dest="mode",
        help="Extract per-checkpoint stats to CSV+JSON",
    )
    mode.add_argument(
        "--snapshot",
        action="store_const",
        const=MODE_SNAPSHOT,
        dest="mode",
        help="Spatial plots for one checkpoint",
    )
    mode.add_argument(
        "--timeseries",
        action="store_const",
        const=MODE_TIMESERIES,
        dest="mode",
        help="Trend plots from a summary CSV",
    )

    parser.add_argument(
        "--checkpoint-dir",
        type=str,
        default=None,
        help="Directory containing step_*.h5 restart files (summarize)",
    )
    parser.add_argument(
        "--checkpoint",
        type=str,
        default=None,
        help="Single checkpoint HDF5 path (snapshot)",
    )
    parser.add_argument(
        "--summary",
        type=str,
        default=None,
        help="Summary CSV from --summarize (timeseries)",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output CSV path (summarize) or directory (snapshot/timeseries)",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Keep every Nth checkpoint when summarizing (default: 1)",
    )
    parser.add_argument(
        "--max-checkpoints",
        type=int,
        default=None,
        help="Optional cap on number of checkpoints to process",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress lines on stderr/stdout",
    )
    return parser


def _select_checkpoints(args: argparse.Namespace) -> list[Path]:
    if not args.checkpoint_dir:
        raise SystemExit("--summarize requires --checkpoint-dir")
    paths = discover_checkpoints(args.checkpoint_dir)
    if args.stride < 1:
        raise SystemExit("--stride must be >= 1")
    paths = paths[:: args.stride]
    if args.max_checkpoints is not None:
        if args.max_checkpoints < 1:
            raise SystemExit("--max-checkpoints must be >= 1")
        paths = paths[: args.max_checkpoints]
    if not paths:
        raise SystemExit(f"no step_*.h5 checkpoints found in {args.checkpoint_dir}")
    return paths


def _meta_path_for_csv(csv_path: Path) -> Path:
    return csv_path.with_name(csv_path.stem + "_meta.json")


def run_summarize(args: argparse.Namespace) -> int:
    paths = _select_checkpoints(args)
    rows, skips = extract_many(paths, progress=not args.quiet)
    if not rows:
        raise SystemExit("no checkpoints could be read")

    frame = pd.DataFrame(rows)
    if "step" in frame.columns:
        frame = frame.sort_values("step").reset_index(drop=True)

    out_csv = prepare_output_file(validate_path_syntax(args.output))
    frame.to_csv(out_csv, index=False)

    meta: dict[str, Any] = {
        "mode": MODE_SUMMARIZE,
        "checkpoint_dir": str(validate_checkpoint_dir(args.checkpoint_dir)),
        "n_requested": len(paths),
        "n_written": len(rows),
        "n_skipped": len(skips),
        "skipped": skips,
        "columns": list(frame.columns),
        "stride": args.stride,
        "csv": str(out_csv),
    }
    meta_path = prepare_output_file(_meta_path_for_csv(out_csv))
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    if not args.quiet:
        print(f"wrote {len(rows)} rows → {out_csv}", flush=True)
        print(f"wrote metadata → {meta_path}", flush=True)
    return 0


def run_snapshot(args: argparse.Namespace) -> int:
    if not args.checkpoint:
        raise SystemExit("--snapshot requires --checkpoint")
    checkpoint = validate_input_path(args.checkpoint)
    out_dir = prepare_output_directory(args.output)
    written = write_snapshot_figures(checkpoint, out_dir)
    if not args.quiet:
        print(f"wrote {len(written)} figures → {out_dir}", flush=True)
        for path in written:
            print(f"  {path.name}", flush=True)
    if not written:
        print("warning: no figures produced (missing agents/grid layers?)", flush=True)
    return 0


def run_timeseries(args: argparse.Namespace) -> int:
    if not args.summary:
        raise SystemExit("--timeseries requires --summary")
    summary_path = validate_input_path(args.summary)
    frame = pd.read_csv(summary_path)
    if frame.empty:
        raise SystemExit(f"summary CSV is empty: {summary_path}")
    out_dir = prepare_output_directory(args.output)
    written = write_timeseries_figures(frame, out_dir)
    if not args.quiet:
        print(f"wrote {len(written)} figures → {out_dir}", flush=True)
        for path in written:
            print(f"  {path.name}", flush=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    handlers = {
        MODE_SUMMARIZE: run_summarize,
        MODE_SNAPSHOT: run_snapshot,
        MODE_TIMESERIES: run_timeseries,
    }
    try:
        return handlers[args.mode](args)
    except PathValidationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except ImportError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except (OSError, ValueError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
