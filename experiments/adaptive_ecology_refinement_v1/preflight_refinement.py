#!/usr/bin/env python3
"""Fail-closed planning/deployment preflight for the adaptive ecology refinement."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
CONTRACT_PATH = ROOT / "refinement_contract.json"
CONTRACT = json.loads(CONTRACT_PATH.read_text())
GENERATED = ROOT / "generated"
EXEC_PLACEHOLDER = CONTRACT["execution_source_sha_policy"]["planning_placeholder"]
IMAGE_PLACEHOLDER = CONTRACT["digest_policy"]["planning_placeholder"]
SHA40 = re.compile(r"[0-9a-f]{40}")
IMAGE = re.compile(CONTRACT["digest_policy"]["image_digest_pattern"])
SEEDS = [int(s) for s in CONTRACT["design"]["seeds"]]
AMPLITUDES = [
    float(a) for a in CONTRACT["design"]["axes"]["bacteriocin.mucin_charge.amplitude"]
]
NULL_AMPLITUDE = float(CONTRACT["design"]["null_amplitude"])
JOBS = int(CONTRACT["design"]["jobs"])
RUNTIME = CONTRACT["runtime"]
PACKAGE_FILES = [
    "refinement_contract.json",
    "refinement_decision_record.json",
    "prepare_refinement.py",
    "preflight_refinement.py",
    "aws_commands_refinement.py",
    "analyze_refinement.py",
    "README.md",
]
ERRORS: list[str] = []
WARNINGS: list[str] = []


def fail(message: str) -> None:
    ERRORS.append(message)


def warn(message: str) -> None:
    WARNINGS.append(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def approx(value, expected, tol=1e-12) -> bool:
    try:
        scale = max(1.0, abs(float(expected)))
        return math.isclose(
            float(value), float(expected), rel_tol=tol, abs_tol=tol * scale
        )
    except (TypeError, ValueError):
        return False


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(REPO), *args], text=True, stderr=subprocess.STDOUT
    ).strip()


def git_head() -> str | None:
    try:
        inside = git("rev-parse", "--is-inside-work-tree")
        return git("rev-parse", "HEAD") if inside == "true" else None
    except (OSError, subprocess.CalledProcessError):
        return None


def check_ancestry(head: str) -> None:
    ancestor = CONTRACT["lineage"]["required_ancestor_sha"]
    try:
        subprocess.check_call(
            ["git", "-C", str(REPO), "cat-file", "-e", ancestor + "^{commit}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        ancestry = subprocess.call(
            ["git", "-C", str(REPO), "merge-base", "--is-ancestor", ancestor, head],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if ancestry != 0:
            fail(f"required ancestor {ancestor} is not an ancestor of HEAD {head}")
    except (OSError, subprocess.CalledProcessError) as exc:
        fail(f"cannot verify required ancestry of {ancestor}: {exc}")


def check_package_tracked_clean() -> None:
    rel = ROOT.relative_to(REPO)
    for name in PACKAGE_FILES:
        try:
            git("ls-files", "--error-unmatch", "--", (rel / name).as_posix())
        except (OSError, subprocess.CalledProcessError):
            fail(
                "refinement package file is not tracked at HEAD: "
                f"{(rel / name).as_posix()}"
            )
    try:
        dirty = git(
            "status", "--porcelain", "--",
            *[(rel / name).as_posix() for name in PACKAGE_FILES],
        )
        if dirty:
            fail("tracked refinement package differs from HEAD:\n" + dirty)
    except (OSError, subprocess.CalledProcessError) as exc:
        fail(f"cannot verify tracked refinement package cleanliness: {exc}")


def validate_git(deployment: bool, supplied_sha: str | None) -> str | None:
    head = git_head()
    if not deployment:
        if head is None:
            warn(
                "planning mode: codeload/non-git tree cannot prove HEAD, ancestry, "
                "or cleanliness; deployment remains blocked"
            )
        else:
            warn(
                "planning mode: git identity is informational only; rerun with "
                "--deployment before command generation"
            )
        return head
    if not supplied_sha or not SHA40.fullmatch(supplied_sha):
        fail("deployment requires --execution-source-sha as full lowercase 40-hex")
    if head is None:
        fail("deployment requires a normal git checkout; codeload is planning-only")
        return head
    if supplied_sha != head:
        fail(f"supplied execution SHA {supplied_sha!r} != git HEAD {head!r}")
    check_ancestry(head)
    check_package_tracked_clean()
    return head


def check_fixed_scalars(rid: str, cfg: dict, amp: float) -> None:
    exact = {
        "total_time": 21600.0,
        "bio_dt": 60.0,
        "output_interval": 300.0,
        "domain_x": 1e-4,
        "domain_y": 1e-4,
        "domain_z": 1e-4,
        "grid_dx": 2e-6,
        "kd_corrinoid_btuB": 1e-4,
        "kd_colicinE_btuB": 5e-7,
        "b12_initial_conc": 1e-3,
        "bacteriocin.mucin_charge.amplitude": amp,
    }
    for key, want in exact.items():
        if not approx(cfg.get(key), want):
            fail(f"{rid}: {key}={cfg.get(key)!r}, expected {want!r}")
    for key, want in {"gpu_enabled": True, "gpu_device_id": 0}.items():
        if cfg.get(key) is not want:
            fail(f"{rid}: {key}={cfg.get(key)!r}, expected {want!r}")
    if cfg.get("fixes") != ["metabolism", "bacteriocin", "receptor", "mechanics"]:
        fail(f"{rid}: fixes must be metabolism/bacteriocin/receptor/mechanics")
    if cfg.get("metabolism.uptake_limit") != "delivery":
        fail(f"{rid}: ecological stream requires metabolism.uptake_limit=delivery")
    if (
        cfg.get("chemistry.toxin_evaluation") != "grid"
        or cfg.get("chemistry.toxin_lumping") != "per_receptor"
    ):
        fail(f"{rid}: chemistry grid/per_receptor settings drift")
    if (
        cfg.get("initial_population.placement") != "z_slab"
        or not approx(cfg.get("initial_population.z_min"), 0.0)
        or not approx(cfg.get("initial_population.z_max"), 1e-4)
    ):
        fail(f"{rid}: required full-depth z_slab placement drift")


def check_strains(rid: str, arm: str, strains: list) -> None:
    if len(strains) != 2 or [s.get("count") for s in strains] != [60, 60]:
        fail(f"{rid}: expected two strains of 60 cells each")
    for s in strains:
        if not approx(s.get("mu_max"), 5.5e-4):
            fail(f"{rid}: every strain must have mu_max=5.5e-4")
        if "receptor_expression" in s:
            fail(f"{rid}: strains must not carry receptor_expression")
    if len(strains) == 2:
        expected_plasmids = ["ColE1"] if arm == "producer" else []
        if strains[0].get("plasmids") != expected_plasmids or strains[1].get(
            "plasmids"
        ) != []:
            fail(f"{rid}: producer/null plasmid identity drift")


def check_schedule_and_meta(
    entry: dict, cfg: dict, expected_sha: str
) -> None:
    rid = entry.get("run_id", "<missing>")
    amp = float(entry.get("amplitude", -1))
    arm = entry.get("arm")
    schedule = (cfg.get("hdf5") or {}).get("schedule") or {}
    if schedule != RUNTIME["hdf5_schedule"]:
        fail(f"{rid}: HDF5 schedule drift: {schedule!r}")
    meta = cfg.get("_refinement") or {}
    if (
        meta.get("refinement_id") != CONTRACT["refinement_id"]
        or meta.get("run_id") != rid
        or meta.get("arm") != arm
        or int(meta.get("seed", -1)) != int(entry.get("seed", -2))
        or not approx(meta.get("amplitude"), amp)
        or meta.get("execution_source_sha") != expected_sha
    ):
        fail(f"{rid}: _refinement metadata does not match manifest/deployment identity")


def validate_config(entry: dict, cfg: dict, expected_sha: str) -> None:
    rid = entry.get("run_id", "<missing>")
    amp = float(entry.get("amplitude", -1))
    check_fixed_scalars(rid, cfg, amp)
    check_strains(rid, entry.get("arm"), cfg.get("initial_strains") or [])
    check_schedule_and_meta(entry, cfg, expected_sha)


def check_run_identity(index: int, entry: dict) -> str:
    rid = str(entry.get("run_id", ""))
    if entry.get("array_index") != index:
        fail(f"array indices are not contiguous at position {index}")
    if (
        not rid.startswith("C_amp")
        or "SS_amp" in rid
        or "single_source" in rid
    ):
        fail(f"{rid!r}: only ecological C_amp run IDs are allowed")
    return rid


def check_input_file(index: int, entry: dict, expected_sha: str) -> None:
    rid = str(entry.get("run_id", ""))
    relraw = entry.get("input_relpath")
    if (
        not isinstance(relraw, str)
        or "\\" in relraw
        or PurePosixPath(relraw).as_posix() != relraw
    ):
        fail(f"{rid}: input_relpath is not portable POSIX syntax: {relraw!r}")
        return
    expected_rel = f"generated/jobs/{index}/input.json"
    if relraw != expected_rel:
        fail(f"{rid}: input_relpath {relraw!r} != {expected_rel!r}")
    path = ROOT / PurePosixPath(relraw)
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        fail(f"{rid}: input path does not resolve on POSIX: {exc}")
        return
    if ROOT.resolve() not in resolved.parents:
        fail(f"{rid}: input path escapes refinement root")
    if digest(path) != entry.get("input_sha256"):
        fail(f"{rid}: input SHA-256 mismatch")
    try:
        cfg = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"{rid}: input JSON unreadable: {exc}")
        return
    validate_config(entry, cfg, expected_sha)


def check_output_and_arm(
    index: int, entry: dict, outputs: set, producers: set, nulls: set
) -> None:
    rid = str(entry.get("run_id", ""))
    out = entry.get("output_relpath")
    expected_out = f"generated/results/{index}/output.h5.gz"
    if out != expected_out or "\\" in str(out):
        fail(f"{rid}: output path is not isolated expected path {expected_out!r}")
    if out in outputs:
        fail(f"{rid}: duplicate output path {out!r}")
    outputs.add(out)
    pair = (float(entry.get("amplitude", -1)), int(entry.get("seed", -1)))
    if entry.get("arm") == "producer":
        producers.add(pair)
    elif entry.get("arm") == "plasmid_free_null":
        nulls.add(pair)
    else:
        fail(f"{rid}: invalid arm {entry.get('arm')!r}")


def validate_runs(manifest: dict, expected_sha: str) -> None:
    runs = manifest.get("runs") or []
    if len(runs) != JOBS:
        fail(f"manifest has {len(runs)} jobs, expected {JOBS}")
    wanted_producers = {(a, s) for a in AMPLITUDES for s in SEEDS}
    wanted_nulls = {(NULL_AMPLITUDE, s) for s in SEEDS}
    seen_producers: set = set()
    seen_nulls: set = set()
    outputs: set = set()
    for index, entry in enumerate(runs):
        check_run_identity(index, entry)
        check_input_file(index, entry, expected_sha)
        check_output_and_arm(index, entry, outputs, seen_producers, seen_nulls)
    if seen_producers != wanted_producers:
        fail(f"producer amplitude/seed set mismatch: {sorted(seen_producers)}")
    if seen_nulls != wanted_nulls:
        fail(f"plasmid-free null set mismatch: {sorted(seen_nulls)}")


def load_manifest() -> dict:
    try:
        return json.loads((GENERATED / "manifest.json").read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"generated deployment manifest missing/unreadable: {exc}")
        return {}


def check_manifest_header(manifest: dict, deployment: bool, expected_sha: str) -> None:
    if manifest.get("schema_version") != 2:
        fail("manifest schema_version must be 2")
    if manifest.get("refinement_id") != CONTRACT["refinement_id"]:
        fail("manifest refinement_id mismatch")
    if manifest.get("refinement_contract_sha256") != digest(CONTRACT_PATH):
        fail("manifest refinement_contract_sha256 mismatch")
    if manifest.get("execution_source_sha") != expected_sha:
        fail(f"manifest execution_source_sha must equal {expected_sha!r}")
    image = manifest.get("container_image_digest")
    if deployment:
        if not isinstance(image, str) or not IMAGE.fullmatch(image):
            fail(
                "deployment manifest requires immutable <repo>@sha256:<64 hex> image"
            )
    elif image != IMAGE_PLACEHOLDER:
        fail("planning manifest must retain the explicit image placeholder")
    else:
        warn("planning mode: image digest is a placeholder; no submission is authorized")


def check_manifest_runtime(manifest: dict) -> None:
    if (
        manifest.get("mpi_ranks") != 1
        or manifest.get("gpu_required") is not True
        or manifest.get("attempt_timeout_s") != 7200
    ):
        fail("manifest must require one MPI rank, one GPU, and 7200 s attempt timeout")
    placement = manifest.get("chemistry_placement")
    if (
        RUNTIME.get("chemistry_placement") != "device_delivery"
        or placement != "device_delivery"
    ):
        fail(
            "ecological refinement must require chemistry_placement=device_delivery; "
            "device, host, and host_forced_delivery are forbidden"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deployment", action="store_true")
    parser.add_argument("--execution-source-sha")
    args = parser.parse_args()
    if args.execution_source_sha and not args.deployment:
        fail("--execution-source-sha is accepted only with --deployment")
    head = validate_git(args.deployment, args.execution_source_sha)

    manifest = load_manifest()
    expected_sha = args.execution_source_sha if args.deployment else EXEC_PLACEHOLDER
    check_manifest_header(manifest, args.deployment, expected_sha)
    check_manifest_runtime(manifest)
    validate_runs(manifest, expected_sha)

    result = {
        "status": "FAIL" if ERRORS else "PASS",
        "mode": "deployment" if args.deployment else "planning",
        "git_head": head,
        "execution_source_sha": manifest.get("execution_source_sha"),
        "runs_checked": len(manifest.get("runs") or []),
        "errors": ERRORS,
        "warnings": WARNINGS,
    }
    print(json.dumps(result, indent=2))
    return 1 if ERRORS else 0


if __name__ == "__main__":
    raise SystemExit(main())
