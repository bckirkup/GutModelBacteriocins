#!/usr/bin/env python3
"""Generate the adaptive ecology refinement inputs. Never runs the model."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO / "python"))
from gut_ibm_tools.transport_metrics import (
    COLE1_LIBRARY_BURST_SIZE,
    COLE1_LIBRARY_DIFF_COEFF,
    COLE1_LIBRARY_PI,
)

CONTRACT_FILENAME = "refinement_contract.json"
CONTRACT = json.loads((ROOT / CONTRACT_FILENAME).read_text())
EXEC_PLACEHOLDER = CONTRACT["execution_source_sha_policy"]["planning_placeholder"]
IMAGE_PLACEHOLDER = CONTRACT["digest_policy"]["planning_placeholder"]
SHA40 = re.compile(r"[0-9a-f]{40}")
IMAGE_RE = re.compile(CONTRACT["digest_policy"]["image_digest_pattern"])
SEEDS = CONTRACT["design"]["seeds"]
AMPLITUDES = CONTRACT["design"]["axes"]["bacteriocin.mucin_charge.amplitude"]
NULL_AMPLITUDE = float(CONTRACT["design"]["null_amplitude"])
JOBS = int(CONTRACT["design"]["jobs"])

DOMAIN_M = 1.0e-4
GRID_DX_M = 2.0e-6
BIO_DT_S = 60.0
TOTAL_TIME_S = 21600.0
STRAIN_COUNT = 60
STRAIN_MU_MAX = 5.5e-4
FIXES = ["metabolism", "bacteriocin", "receptor", "mechanics"]

PACKAGE_FILES = [
    CONTRACT_FILENAME,
    "refinement_decision_record.json",
    "prepare_refinement.py",
    "preflight_refinement.py",
    "aws_commands_refinement.py",
    "analyze_refinement.py",
    "README.md",
]


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(REPO), *args], text=True, stderr=subprocess.STDOUT
    ).strip()


def dump(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(obj, indent=2, sort_keys=True) + "\n").encode("utf-8")
    handle, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_deployment_identity(execution_sha: str, image_digest: str) -> None:
    if not SHA40.fullmatch(execution_sha or ""):
        raise SystemExit(
            "REFUSED: --deployment requires --execution-source-sha as a full "
            "lowercase 40-hex commit"
        )
    if not IMAGE_RE.fullmatch(image_digest or ""):
        raise SystemExit(
            "REFUSED: --deployment requires --image-digest <repo>@sha256:<64 hex>"
        )
    try:
        inside = git("rev-parse", "--is-inside-work-tree")
        head = git("rev-parse", "HEAD")
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(
            f"REFUSED: deployment generation requires a git checkout: {exc}"
        ) from exc
    if inside != "true" or head != execution_sha:
        raise SystemExit(
            f"REFUSED: --execution-source-sha {execution_sha} does not match git HEAD {head}"
        )
    ancestor = CONTRACT["lineage"]["required_ancestor_sha"]
    try:
        git("cat-file", "-e", ancestor + "^{commit}")
        ancestry = subprocess.call(
            ["git", "-C", str(REPO), "merge-base", "--is-ancestor", ancestor, head],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(
            f"REFUSED: cannot verify required ancestry of {ancestor}: {exc}"
        ) from exc
    if ancestry != 0:
        raise SystemExit(
            f"REFUSED: required ancestor {ancestor} is not an ancestor of HEAD {head}"
        )
    relative = ROOT.relative_to(REPO)
    dirty = git(
        "status", "--porcelain", "--", *[str(relative / name) for name in PACKAGE_FILES]
    )
    if dirty:
        raise SystemExit(
            "REFUSED: refinement package files differ from HEAD; commit them before "
            "deployment generation:\n" + dirty
        )


def strain(strain_type: int, *, plasmids: list[str]) -> dict:
    return {
        "type": strain_type,
        "count": STRAIN_COUNT,
        "mu_max": STRAIN_MU_MAX,
        "plasmids": plasmids,
        "conjugative": False,
    }


def config(
    *, run_id: str, arm: str, seed: int, amplitude: float, execution_sha: str
) -> dict:
    producer_plasmids = ["ColE1"] if arm == "producer" else []
    return {
        "_comment": [
            (
                "Adaptive ecology refinement: same-image "
                "producer/plasmid-free-null ecological contrast for Gate C."
            ),
            f"Execution source: {execution_sha}.",
            "Runtime /run_provenance/git_sha must match execution_source_sha.",
        ],
        "_refinement": {
            "refinement_id": CONTRACT["refinement_id"],
            "run_id": run_id,
            "arm": arm,
            "seed": seed,
            "amplitude": amplitude,
            "execution_source_sha": execution_sha,
            "cole1_identity": {
                "pI": COLE1_LIBRARY_PI,
                "diff_coeff": COLE1_LIBRARY_DIFF_COEFF,
                "burst_size": COLE1_LIBRARY_BURST_SIZE,
            },
        },
        "total_time": TOTAL_TIME_S,
        "bio_dt": BIO_DT_S,
        "output_interval": 300.0,
        "seed": seed,
        "domain_x": DOMAIN_M,
        "domain_y": DOMAIN_M,
        "domain_z": DOMAIN_M,
        "grid_dx": GRID_DX_M,
        "mucus_thickness": DOMAIN_M,
        "radial_turnover": 1.0e9,
        "distal_transit": 1.0e9,
        "peristaltic_enabled": False,
        "crypts_enabled": False,
        "motility.enabled": False,
        "carbon_z_gradient": False,
        "carbon.boundary_conc": 0.05,
        "metabolism.uptake_limit": "delivery",
        "oxygen.k_ROS": 0.0,
        "dysbiosis_threshold": 1.0e10,
        "gpu_enabled": True,
        "gpu_device_id": 0,
        "chemistry.toxin_evaluation": "grid",
        "chemistry.toxin_lumping": "per_receptor",
        "initial_population.placement": "z_slab",
        "initial_population.z_min": 0.0,
        "initial_population.z_max": DOMAIN_M,
        "fixes": list(FIXES),
        "kd_colicinE_btuB": 5.0e-7,
        "kd_corrinoid_btuB": 1.0e-4,
        "b12_initial_conc": 1.0e-3,
        "bacteriocin.mucin_charge.amplitude": amplitude,
        "initial_strains": [
            strain(1, plasmids=producer_plasmids),
            strain(2, plasmids=[]),
        ],
        "hdf5_file": "output.h5",
        "hdf5": {
            "enabled": True,
            "compression": "gzip",
            "compression_level": 4,
            "schedule": dict(CONTRACT["runtime"]["hdf5_schedule"]),
        },
    }


def planned_runs(execution_sha: str) -> list[tuple[str, dict]]:
    runs: list[tuple[str, dict]] = []
    for amplitude in AMPLITUDES:
        for seed in SEEDS:
            run_id = f"C_amp{int(amplitude)}_producer_s{seed}"
            runs.append(
                (
                    run_id,
                    config(
                        run_id=run_id,
                        arm="producer",
                        seed=seed,
                        amplitude=float(amplitude),
                        execution_sha=execution_sha,
                    ),
                )
            )
    for seed in SEEDS:
        run_id = f"C_amp{int(NULL_AMPLITUDE)}_plasmid_free_null_s{seed}"
        runs.append(
            (
                run_id,
                config(
                    run_id=run_id,
                    arm="plasmid_free_null",
                    seed=seed,
                    amplitude=NULL_AMPLITUDE,
                    execution_sha=execution_sha,
                ),
            )
        )
    return runs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deployment", action="store_true")
    parser.add_argument("--execution-source-sha")
    parser.add_argument("--image-digest", default=IMAGE_PLACEHOLDER)
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()

    execution_sha = args.execution_source_sha or EXEC_PLACEHOLDER
    if args.deployment:
        validate_deployment_identity(execution_sha, args.image_digest)
    elif args.execution_source_sha or args.image_digest != IMAGE_PLACEHOLDER:
        raise SystemExit(
            "REFUSED: a real execution SHA or image digest requires --deployment; "
            "planning generation uses explicit placeholders"
        )

    generated = ROOT / "generated"
    if args.clean and generated.exists():
        shutil.rmtree(generated)

    entries = []
    for index, (run_id, cfg) in enumerate(planned_runs(execution_sha)):
        input_path = generated / "jobs" / str(index) / "input.json"
        dump(input_path, cfg)
        entries.append(
            {
                "array_index": index,
                "run_id": run_id,
                "arm": cfg["_refinement"]["arm"],
                "seed": cfg["_refinement"]["seed"],
                "amplitude": cfg["_refinement"]["amplitude"],
                "input_relpath": input_path.relative_to(ROOT).as_posix(),
                "input_sha256": digest(input_path),
                "output_relpath": f"generated/results/{index}/output.h5.gz",
            }
        )

    if len(entries) != JOBS:
        raise SystemExit(
            f"REFUSED: generated {len(entries)} runs but the contract declares {JOBS}"
        )

    dump(
        generated / "manifest.json",
        {
            "schema_version": 2,
            "refinement_id": CONTRACT["refinement_id"],
            "refinement_contract_sha256": digest(ROOT / CONTRACT_FILENAME),
            "execution_source_sha": execution_sha,
            "container_image_digest": args.image_digest,
            "mpi_ranks": CONTRACT["runtime"]["mpi_ranks"],
            "gpu_required": CONTRACT["runtime"]["gpu_required"],
            "chemistry_placement": CONTRACT["runtime"]["chemistry_placement"],
            "attempt_timeout_s": CONTRACT["runtime"]["attempt_timeout_s"],
            "gate": CONTRACT["gate"]["id"],
            "runs": entries,
        },
    )
    print(
        json.dumps(
            {
                "generated_runs": len(entries),
                "execution_source_sha": execution_sha,
                "container_image_digest": args.image_digest,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
