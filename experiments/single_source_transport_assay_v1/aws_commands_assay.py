#!/usr/bin/env python3
"""Validate deployment identity, then print (never execute) assay AWS commands."""
from __future__ import annotations
import argparse, json, re, shlex, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
SHA40 = re.compile(r"[0-9a-f]{40}")
IMAGE = re.compile(r"[^\s@]+@sha256:[0-9a-f]{64}")
S3 = re.compile(r"s3://[^\s/]+/.+")

def refuse(message: str) -> None: raise SystemExit("REFUSED: " + message)
def git(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(REPO), *args], text=True, stderr=subprocess.STDOUT).strip()

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--execution-source-sha", required=True)
    ap.add_argument("--image-uri", required=True, help="immutable image URI <repo>@sha256:<64 hex>")
    ap.add_argument("--input-prefix", required=True, help="isolated s3://.../single-source-assay.../inputs prefix")
    ap.add_argument("--output-prefix", required=True, help="isolated s3://.../single-source-assay.../outputs prefix")
    ap.add_argument("--job-name", default="gutibm-single-source-assay-v1-adaptive")
    ap.add_argument("--job-queue", required=True)
    ap.add_argument("--job-definition", required=True)
    ap.add_argument("--authorization-file", type=Path, help="optional operator authorization JSON, validated if supplied")
    args = ap.parse_args()
    if not SHA40.fullmatch(args.execution_source_sha): refuse("execution SHA must be full lowercase 40-hex")
    if not IMAGE.fullmatch(args.image_uri): refuse("image URI must be immutable <repo>@sha256:<64 hex>")
    for label, value in (("input", args.input_prefix), ("output", args.output_prefix)):
        lowered = value.lower()
        if not S3.fullmatch(value): refuse(f"{label} prefix must be a non-root s3:// URI")
        if any(token in lowered for token in ("receptor-v1/stage", "stage_c", "stage-c", "/ecological", "c_amp")):
            refuse(f"{label} prefix looks like a stale ecological campaign path")
        if "single-source" not in lowered or "assay" not in lowered: refuse(f"{label} prefix must identify an isolated single-source assay path")
    if args.input_prefix.rstrip("/") == args.output_prefix.rstrip("/"): refuse("input and output prefixes must be isolated")
    if not args.job_name.startswith("gutibm-single-source-assay") or "c-amp" in args.job_name.lower(): refuse("job name must identify only the single-source assay")
    try: head = git("rev-parse", "HEAD")
    except Exception as exc: refuse(f"normal git checkout required: {exc}")
    if head != args.execution_source_sha: refuse(f"execution SHA differs from git HEAD {head}")
    manifest_path = ROOT / "generated" / "manifest.json"
    try: manifest = json.loads(manifest_path.read_text())
    except Exception as exc: refuse(f"deployment manifest missing/unreadable: {exc}")
    if manifest.get("execution_source_sha") != head: refuse("deployment manifest SHA differs from exact HEAD; regenerate it")
    if manifest.get("container_image_digest") != args.image_uri: refuse("--image-uri is not byte-for-byte equal to deployment manifest image")
    if len(manifest.get("runs") or []) != 18: refuse("deployment manifest is not the 18-job adaptive intrinsic assay")
    ids = [str(run.get("run_id", "")) for run in manifest["runs"]]
    if any(not rid.startswith("SS_amp") or "C_amp" in rid for rid in ids): refuse("manifest contains stale ecological/non-assay run IDs")
    # Reuse the complete deployment preflight; capture it so this program emits commands only on PASS.
    check = subprocess.run([sys.executable, str(ROOT / "preflight_assay.py"), "--deployment", "--execution-source-sha", head], text=True, capture_output=True, check=False)
    if check.returncode != 0: refuse("deployment preflight failed:\n" + check.stdout + check.stderr)
    if args.authorization_file:
        try: auth = json.loads(args.authorization_file.read_text())
        except Exception as exc: refuse(f"authorization file unreadable: {exc}")
        if auth.get("authorize_single_source_assay_submission") is not True or auth.get("execution_source_sha") != head or auth.get("container_image_digest") != args.image_uri or auth.get("jobs") != 18:
            refuse("authorization must approve this exact 18-job assay SHA and image")
    q = shlex.quote
    local_jobs = ROOT / "generated" / "jobs"
    in_prefix = args.input_prefix.rstrip("/")
    out_prefix = args.output_prefix.rstrip("/")
    env = [
        {"name":"INPUT_S3_PREFIX","value":in_prefix},
        {"name":"OUTPUT_S3_PREFIX","value":out_prefix},
        {"name":"MPI_RANKS","value":"1"},
        {"name":"REQUIRE_GPU","value":"1"},
        {"name":"GPU_DEVICE_ID","value":"0"},
    ]
    submit = ["aws","batch","submit-job","--job-name",args.job_name,"--job-queue",args.job_queue,"--job-definition",args.job_definition,"--array-properties","size=18","--timeout","attemptDurationSeconds=3600","--container-overrides",json.dumps({"environment":env},separators=(",",":"))]
    print("# REVIEW ONLY. This program did not execute AWS commands.")
    print(f"# execution_source_sha={head}")
    print(f"# image_uri={args.image_uri}")
    print(f"# Verify job definition {args.job_definition} resolves exactly to that image URI.")
    print(f"aws s3 cp {q(str(local_jobs))} {q(in_prefix + '/')} --recursive --exclude '*' --include '*/input.json'")
    print(" ".join(q(part) for part in submit))
    return 0

if __name__ == "__main__": sys.exit(main())
