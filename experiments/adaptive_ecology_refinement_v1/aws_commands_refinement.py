#!/usr/bin/env python3
"""Validate deployment identity, then print (never execute) refinement AWS commands."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO / "python"))
from gut_ibm_tools.path_utils import PathValidationError, validate_input_path

CONTRACT = json.loads((ROOT / "refinement_contract.json").read_text())
JOBS = int(CONTRACT["design"]["jobs"])
SHA40 = re.compile(r"[0-9a-f]{40}")
IMAGE = re.compile(CONTRACT["digest_policy"]["image_digest_pattern"])
S3 = re.compile(r"s3://[^\s/]+/.+")
STALE_TOKENS = ("single-source", "ss_amp", "receptor-v1/stage", "stage_c", "stage-c")
PREFIX_TOKEN = "adaptive-ecology-refinement"
JOB_NAME_PREFIX = "gutibm-adaptive-ecology-refinement"


def refuse(message: str) -> None:
    raise SystemExit("REFUSED: " + message)


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(REPO), *args], text=True, stderr=subprocess.STDOUT
    ).strip()


def load_manifest() -> dict:
    try:
        return json.loads((ROOT / "generated" / "manifest.json").read_text())
    except (OSError, json.JSONDecodeError) as exc:
        refuse(f"deployment manifest missing/unreadable: {exc}")


def validate_prefix(label: str, value: str) -> None:
    lowered = value.lower()
    if not S3.fullmatch(value):
        refuse(f"{label} prefix must be a non-root s3:// URI")
    if PREFIX_TOKEN not in lowered:
        refuse(
            f"{label} prefix must identify an isolated {PREFIX_TOKEN} path"
        )
    if any(token in lowered for token in STALE_TOKENS):
        refuse(f"{label} prefix looks like an intrinsic assay or stale stage-C path")


def _validate_args(args) -> None:
    if not SHA40.fullmatch(args.execution_source_sha):
        refuse("execution SHA must be full lowercase 40-hex")
    if not IMAGE.fullmatch(args.image_uri):
        refuse("image URI must be immutable <repo>@sha256:<64 hex>")
    validate_prefix("input", args.input_prefix)
    validate_prefix("output", args.output_prefix)
    if args.input_prefix.rstrip("/") == args.output_prefix.rstrip("/"):
        refuse("input and output prefixes must be isolated")
    job_name = args.job_name.lower()
    if not job_name.startswith(JOB_NAME_PREFIX) or any(
        token in job_name for token in ("single-source", "ss-amp")
    ):
        refuse("job name must identify only the adaptive ecology refinement")


def _validate_against_deployment(args, manifest: dict, head: str) -> None:
    if manifest.get("execution_source_sha") != head:
        refuse("deployment manifest SHA differs from exact HEAD; regenerate it")
    if manifest.get("container_image_digest") != args.image_uri:
        refuse("--image-uri is not byte-for-byte equal to deployment manifest image")
    runs = manifest.get("runs") or []
    if len(runs) != JOBS:
        refuse(f"deployment manifest is not the {JOBS}-job ecological refinement")
    ids = [str(run.get("run_id", "")) for run in runs]
    if any(not rid.startswith("C_amp") or "SS_amp" in rid for rid in ids):
        refuse("manifest contains non-ecological run IDs")


def _run_preflight(head: str) -> None:
    # Reuse the complete deployment preflight; capture it so this program emits
    # commands only on PASS.
    check = subprocess.run(
        [
            sys.executable,
            str(ROOT / "preflight_refinement.py"),
            "--deployment",
            "--execution-source-sha",
            head,
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    if check.returncode != 0:
        refuse("deployment preflight failed:\n" + check.stdout + check.stderr)


def _load_authorization(args, image_uri: str, head: str) -> None:
    if not args.authorization_file:
        return
    try:
        auth = json.loads(validate_input_path(args.authorization_file).read_text())
    except (OSError, json.JSONDecodeError, PathValidationError) as exc:
        refuse(f"authorization file unreadable: {exc}")
    if (
        auth.get("authorize_adaptive_ecology_refinement_submission") is not True
        or auth.get("execution_source_sha") != head
        or auth.get("container_image_digest") != image_uri
        or auth.get("jobs") != JOBS
    ):
        refuse(
            f"authorization must approve this exact {JOBS}-job refinement SHA "
            "and image"
        )


def _print_commands(args, head: str) -> None:
    q = shlex.quote
    local_jobs = ROOT / "generated" / "jobs"
    in_prefix = args.input_prefix.rstrip("/")
    out_prefix = args.output_prefix.rstrip("/")
    env = [
        {"name": "INPUT_S3_PREFIX", "value": in_prefix},
        {"name": "OUTPUT_S3_PREFIX", "value": out_prefix},
        {"name": "MPI_RANKS", "value": "1"},
        {"name": "REQUIRE_GPU", "value": "1"},
        {"name": "GPU_DEVICE_ID", "value": "0"},
    ]
    submit = [
        "aws", "batch", "submit-job",
        "--job-name", args.job_name,
        "--job-queue", args.job_queue,
        "--job-definition", args.job_definition,
        "--array-properties", f"size={JOBS}",
        "--timeout", "attemptDurationSeconds=7200",
        "--container-overrides",
        json.dumps({"environment": env}, separators=(",", ":")),
    ]
    print("# REVIEW ONLY. This program did not execute AWS commands.")
    print(f"# execution_source_sha={head}")
    print(f"# image_uri={args.image_uri}")
    print(
        f"# Verify job definition {args.job_definition} resolves exactly to that "
        "image URI."
    )
    print(
        f"aws s3 cp {q(str(local_jobs))} {q(in_prefix + '/')} --recursive "
        "--exclude '*' --include '*/input.json'"
    )
    print(" ".join(q(part) for part in submit))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--execution-source-sha", required=True)
    ap.add_argument(
        "--image-uri", required=True, help="immutable image URI <repo>@sha256:<64 hex>"
    )
    ap.add_argument(
        "--input-prefix",
        required=True,
        help="isolated s3://.../adaptive-ecology-refinement.../inputs prefix",
    )
    ap.add_argument(
        "--output-prefix",
        required=True,
        help="isolated s3://.../adaptive-ecology-refinement.../outputs prefix",
    )
    ap.add_argument("--job-name", default="gutibm-adaptive-ecology-refinement-v1")
    ap.add_argument("--job-queue", required=True)
    ap.add_argument("--job-definition", required=True)
    ap.add_argument(
        "--authorization-file",
        type=Path,
        help="optional operator authorization JSON, validated if supplied",
    )
    args = ap.parse_args()

    _validate_args(args)
    try:
        head = git("rev-parse", "HEAD")
    except (OSError, subprocess.CalledProcessError) as exc:
        refuse(f"normal git checkout required: {exc}")
    if head != args.execution_source_sha:
        refuse(f"execution SHA differs from git HEAD {head}")
    manifest = load_manifest()
    _validate_against_deployment(args, manifest, head)
    _run_preflight(head)
    _load_authorization(args, args.image_uri, head)
    _print_commands(args, head)
    return 0


if __name__ == "__main__":
    sys.exit(main())
