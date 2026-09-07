#!/usr/bin/env python3
"""Fail-closed planning/deployment preflight for the intrinsic transport assay."""
from __future__ import annotations
import argparse, hashlib, json, math, re, subprocess, sys
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
CONTRACT_PATH = ROOT / "assay_contract.json"
CONTRACT = json.loads(CONTRACT_PATH.read_text())
GENERATED = ROOT / "generated"
EXEC_PLACEHOLDER = "EXECUTION_SOURCE_SHA_REQUIRED_AFTER_ASSAY_COMMIT"
IMAGE_PLACEHOLDER = "REQUIRED_BEFORE_SUBMISSION"
SHA40 = re.compile(r"[0-9a-f]{40}")
IMAGE = re.compile(r"[^\s@]+@sha256:[0-9a-f]{64}")
PACKAGE_FILES = [
    "assay_contract.json", "assay_decision_record.json", "prepare_assay.py",
    "preflight_assay.py", "aws_commands_assay.py", "analyze_assay.py", "README.md",
]
ERRORS: list[str] = []
WARNINGS: list[str] = []

def fail(message: str) -> None: ERRORS.append(message)
def warn(message: str) -> None: WARNINGS.append(message)
def digest(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def approx(value, expected, tol=1e-12) -> bool:
    try: return math.isclose(float(value), float(expected), rel_tol=tol, abs_tol=tol * max(1.0, abs(float(expected))))
    except (TypeError, ValueError): return False

def git(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(REPO), *args], text=True, stderr=subprocess.STDOUT).strip()

def git_head() -> str | None:
    try:
        return git("rev-parse", "HEAD") if git("rev-parse", "--is-inside-work-tree") == "true" else None
    except Exception:
        return None

def validate_git(deployment: bool, supplied_sha: str | None) -> str | None:
    head = git_head()
    ancestor = CONTRACT["lineage"]["required_ancestor_sha"]
    if not deployment:
        if head is None:
            warn("planning mode: codeload/non-git tree cannot prove HEAD, ancestry, or cleanliness; deployment remains blocked")
        else:
            warn("planning mode: git identity is informational only; rerun with --deployment before command generation")
        return head
    if not supplied_sha or not SHA40.fullmatch(supplied_sha):
        fail("deployment requires --execution-source-sha as full lowercase 40-hex")
    if head is None:
        fail("deployment requires a normal git checkout; codeload is planning-only")
        return head
    if supplied_sha != head:
        fail(f"supplied execution SHA {supplied_sha!r} != git HEAD {head!r}")
    try:
        subprocess.check_call(["git", "-C", str(REPO), "cat-file", "-e", ancestor + "^{commit}"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if subprocess.call(["git", "-C", str(REPO), "merge-base", "--is-ancestor", ancestor, head], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            fail(f"required PR417 event-time ancestor {ancestor} is not an ancestor of HEAD {head}")
    except Exception as exc:
        fail(f"cannot verify required PR417 event-time ancestry: {exc}")
    rel = ROOT.relative_to(REPO)
    for name in PACKAGE_FILES:
        try: git("ls-files", "--error-unmatch", "--", (rel / name).as_posix())
        except Exception: fail(f"assay package file is not tracked at HEAD: {(rel/name).as_posix()}")
    try:
        dirty = git("status", "--porcelain", "--", *[(rel / name).as_posix() for name in PACKAGE_FILES])
        if dirty: fail("tracked assay package differs from HEAD:\n" + dirty)
    except Exception as exc: fail(f"cannot verify tracked assay package cleanliness: {exc}")
    return head

def validate_config(entry: dict, cfg: dict, expected_sha: str) -> None:
    rid = entry.get("run_id", "<missing>")
    amp = float(entry.get("amplitude", -1))
    arm = entry.get("arm")
    exact = {
        "total_time": 3600.0, "bio_dt": 60.0, "domain_x": 1e-4,
        "domain_y": 1e-4, "domain_z": 1e-4, "grid_dx": 2e-6,
        "kd_corrinoid_btuB": 1e-4, "b12_initial_conc": 1e-3,
        "burst_release_tau": 300.0, "bacteriocin.mucin_charge.amplitude": amp,
    }
    for key, want in exact.items():
        if not approx(cfg.get(key), want): fail(f"{rid}: {key}={cfg.get(key)!r}, expected {want!r}")
    for key, want in {"gpu_enabled": True, "gpu_device_id": 0}.items():
        if cfg.get(key) is not want: fail(f"{rid}: {key}={cfg.get(key)!r}, expected {want!r}")
    if cfg.get("fixes") != ["bacteriocin", "receptor"]: fail(f"{rid}: fixes must enable only bacteriocin and receptor")
    if any(key == "metabolism.uptake_limit" or key.startswith("metabolism.") for key in cfg):
        fail(f"{rid}: intrinsic assay must not configure metabolism")
    if cfg.get("chemistry.toxin_evaluation") != "grid" or cfg.get("chemistry.toxin_lumping") != "per_receptor":
        fail(f"{rid}: chemistry grid/per_receptor settings drift")
    if cfg.get("initial_population.placement") != "z_slab" or not approx(cfg.get("initial_population.z_min"), 49e-6) or not approx(cfg.get("initial_population.z_max"), 51e-6):
        fail(f"{rid}: required 2 um mid-depth placement band drift")
    strains = cfg.get("initial_strains") or []
    if len(strains) != 2 or [x.get("count") for x in strains] != [1, 3]: fail(f"{rid}: expected one source cell plus three bystanders")
    for strain in strains:
        if not approx(strain.get("mu_max"), 0.0) or strain.get("receptor_expression") != {"BtuB": 0.0}:
            fail(f"{rid}: every cell type must have mu_max=0 and BtuB=0")
    if len(strains) == 2:
        expected_plasmids = ["ColE1"] if arm == "producer" else []
        if strains[0].get("plasmids") != expected_plasmids or strains[1].get("plasmids") != []:
            fail(f"{rid}: producer/null plasmid identity drift")
    schedule = (cfg.get("hdf5") or {}).get("schedule") or {}
    expected_schedule = {"summary":1,"agents":10,"grid":2,"lineage":0,"genome":10,"provenance":1,"grid_species":["bacteriocin_BtuB"]}
    if schedule != expected_schedule: fail(f"{rid}: HDF5 schedule drift: {schedule!r}")
    meta = cfg.get("_assay") or {}
    if meta.get("run_id") != rid or meta.get("arm") != arm or int(meta.get("seed", -1)) != int(entry.get("seed", -2)) or not approx(meta.get("amplitude"), amp) or meta.get("execution_source_sha") != expected_sha:
        fail(f"{rid}: _assay metadata does not match manifest/deployment identity")

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--deployment", action="store_true")
    ap.add_argument("--execution-source-sha")
    args = ap.parse_args()
    if args.execution_source_sha and not args.deployment: fail("--execution-source-sha is accepted only with --deployment")
    head = validate_git(args.deployment, args.execution_source_sha)
    manifest_path = GENERATED / "manifest.json"
    try: manifest = json.loads(manifest_path.read_text())
    except Exception as exc:
        fail(f"generated deployment manifest missing/unreadable: {exc}"); manifest = {}
    expected_sha = args.execution_source_sha if args.deployment else EXEC_PLACEHOLDER
    if manifest.get("schema_version") != 2: fail("manifest schema_version must be 2")
    if manifest.get("assay_id") != CONTRACT["assay_id"]: fail("manifest assay_id mismatch")
    if manifest.get("assay_contract_sha256") != digest(CONTRACT_PATH): fail("manifest assay_contract_sha256 mismatch")
    if manifest.get("execution_source_sha") != expected_sha: fail(f"manifest execution_source_sha must equal {expected_sha!r}")
    image = manifest.get("container_image_digest")
    if args.deployment:
        if not isinstance(image, str) or not IMAGE.fullmatch(image): fail("deployment manifest requires immutable <repo>@sha256:<64 hex> image")
    elif image != IMAGE_PLACEHOLDER:
        fail("planning manifest must retain the explicit image placeholder")
    else: warn("planning mode: image digest is a placeholder; no submission is authorized")
    if manifest.get("mpi_ranks") != 1 or manifest.get("gpu_required") is not True or manifest.get("attempt_timeout_s") != 3600:
        fail("manifest must require one MPI rank, one GPU, and 3600 s attempt timeout")
    if CONTRACT["runtime"].get("chemistry_placement") != "device" or manifest.get("chemistry_placement") != "device":
        fail("intrinsic assay must require physical chemistry_placement=device; host and host_forced_delivery are forbidden")
    runs = manifest.get("runs") or []
    if len(runs) != 18: fail(f"manifest has {len(runs)} jobs, expected 18")
    wanted_producers = {(float(a), int(seed)) for a in [0,15,20,30,60] for seed in CONTRACT["design"]["seeds"]}
    seen_producers, seen_nulls, outputs = set(), set(), set()
    for index, entry in enumerate(runs):
        rid = str(entry.get("run_id", ""))
        if entry.get("array_index") != index: fail(f"array indices are not contiguous at position {index}")
        if not rid.startswith("SS_amp") or "C_amp" in rid: fail(f"{rid!r}: only intrinsic SS_amp IDs are allowed")
        relraw = entry.get("input_relpath")
        if not isinstance(relraw, str) or "\\" in relraw or PurePosixPath(relraw).as_posix() != relraw:
            fail(f"{rid}: input_relpath is not portable POSIX syntax: {relraw!r}"); continue
        expected_rel = f"generated/jobs/{index}/input.json"
        if relraw != expected_rel: fail(f"{rid}: input_relpath {relraw!r} != {expected_rel!r}")
        path = ROOT / PurePosixPath(relraw)
        try: resolved = path.resolve(strict=True)
        except Exception as exc: fail(f"{rid}: input path does not resolve on POSIX: {exc}"); continue
        if ROOT.resolve() not in resolved.parents: fail(f"{rid}: input path escapes assay root")
        if digest(path) != entry.get("input_sha256"): fail(f"{rid}: input SHA-256 mismatch")
        try: cfg = json.loads(path.read_text())
        except Exception as exc: fail(f"{rid}: input JSON unreadable: {exc}"); continue
        validate_config(entry, cfg, expected_sha)
        out = entry.get("output_relpath")
        expected_out = f"generated/results/{index}/output.h5.gz"
        if out != expected_out or "\\" in str(out): fail(f"{rid}: output path is not isolated expected path {expected_out!r}")
        if out in outputs: fail(f"{rid}: duplicate output path {out!r}")
        outputs.add(out)
        pair = (float(entry.get("amplitude", -1)), int(entry.get("seed", -1)))
        if entry.get("arm") == "producer": seen_producers.add(pair)
        elif entry.get("arm") == "toxin_free_null": seen_nulls.add(pair)
        else: fail(f"{rid}: invalid arm {entry.get('arm')!r}")
    if seen_producers != wanted_producers: fail(f"producer amplitude/seed set mismatch: {sorted(seen_producers)}")
    wanted_nulls = {(60.0, int(seed)) for seed in CONTRACT["design"]["seeds"]}
    if seen_nulls != wanted_nulls: fail(f"null set mismatch: {sorted(seen_nulls)}")
    result = {"status":"FAIL" if ERRORS else "PASS", "mode":"deployment" if args.deployment else "planning", "git_head":head, "execution_source_sha":manifest.get("execution_source_sha"), "runs_checked":len(runs), "errors":ERRORS, "warnings":WARNINGS}
    print(json.dumps(result, indent=2))
    return 1 if ERRORS else 0

if __name__ == "__main__": sys.exit(main())
