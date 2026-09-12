#!/usr/bin/env python3
"""Generate the single-source transport assay inputs. Never runs the model."""

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

ASSAY_CONTRACT_FILE = "assay_contract.json"
CONTRACT = json.loads((ROOT / ASSAY_CONTRACT_FILE).read_text())
EXEC_PLACEHOLDER = "EXECUTION_SOURCE_SHA_REQUIRED_AFTER_ASSAY_COMMIT"
IMAGE_PLACEHOLDER = "REQUIRED_BEFORE_SUBMISSION"
SHA40 = re.compile(r"[0-9a-f]{40}")
IMAGE_RE = re.compile(r"[^\s@]+@sha256:[0-9a-f]{64}")
SEEDS = CONTRACT["design"]["seeds"]
AMPLITUDES = CONTRACT["design"]["axes"]["bacteriocin.mucin_charge.amplitude"]

DOMAIN_M = 1.0e-4
GRID_DX_M = 2.0e-6
BIO_DT_S = 60.0
TOTAL_TIME_S = 3600.0
GRID_EVERY_STEPS = 2  # 120 s grid snapshots: dense pairing inside the release window
SOS_BASAL_RATE = 1.0  # the single producer is SOS-induced on its first step
BYSTANDER_COUNT = 3  # keeps the global agent count above the population-stop threshold
PLACEMENT_HALF_BAND_M = 1.0e-6


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
    tracked = [
        ASSAY_CONTRACT_FILE,
        "assay_decision_record.json",
        "prepare_assay.py",
        "preflight_assay.py",
        "aws_commands_assay.py",
        "analyze_assay.py",
        "README.md",
    ]
    relative = ROOT.relative_to(REPO)
    dirty = git(
        "status", "--porcelain", "--", *[str(relative / name) for name in tracked]
    )
    if dirty:
        raise SystemExit(
            "REFUSED: assay package files differ from HEAD; commit them before "
            "deployment generation:\n" + dirty
        )


def strain(strain_type: int, *, plasmids: list[str], count: int) -> dict:
    return {
        "type": strain_type,
        "count": count,
        "mu_max": 0.0,
        "plasmids": plasmids,
        "conjugative": False,
        "receptor_expression": {"BtuB": 0.0},
    }


def config(
    *, run_id: str, arm: str, seed: int, amplitude: float, execution_sha: str
) -> dict:
    producer_plasmids = ["ColE1"] if arm == "producer" else []
    return {
        "_comment": [
            "Single-source transport assay: one identical fixed release source per seed.",
            f"Execution source: {execution_sha}.",
            "Runtime /run_provenance/git_sha must match execution_source_sha.",
        ],
        "_assay": {
            "assay_id": CONTRACT["assay_id"],
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
        "oxygen.k_ROS": 0.0,
        "dysbiosis_threshold": 1.0e10,
        "gpu_enabled": True,
        "gpu_device_id": 0,
        "chemistry.toxin_evaluation": "grid",
        "chemistry.toxin_lumping": "per_receptor",
        "initial_population.placement": "z_slab",
        "initial_population.z_min": 0.5 * DOMAIN_M - PLACEMENT_HALF_BAND_M,
        "initial_population.z_max": 0.5 * DOMAIN_M + PLACEMENT_HALF_BAND_M,
        "fixes": ["bacteriocin", "receptor"],
        "kd_corrinoid_btuB": 1.0e-4,
        "kd_colicinE_btuB": 5.0e-7,
        "b12_initial_conc": 1.0e-3,
        "burst_release_tau": 300.0,
        "bacteriocin.mucin_charge.amplitude": amplitude,
        "sos_basal_rate": SOS_BASAL_RATE,
        "sos_lysis_prob": 0.0,
        "initial_strains": [
            strain(1, plasmids=producer_plasmids, count=1),
            strain(2, plasmids=[], count=BYSTANDER_COUNT),
        ],
        "hdf5_file": "output.h5",
        "hdf5": {
            "enabled": True,
            "compression": "gzip",
            "compression_level": 4,
            "schedule": {
                "summary": 1,
                "agents": 10,
                "grid": GRID_EVERY_STEPS,
                "lineage": 0,
                "genome": 10,
                "provenance": 1,
                "grid_species": ["bacteriocin_BtuB"],
            },
        },
    }


def planned_runs(execution_sha: str) -> list[tuple[str, dict]]:
    runs: list[tuple[str, dict]] = []
    for amplitude in AMPLITUDES:
        for seed in SEEDS:
            run_id = f"SS_amp{amplitude}_producer_s{seed}"
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
    null_amplitude = float(AMPLITUDES[-1])
    for seed in SEEDS:
        run_id = f"SS_amp{int(null_amplitude)}_toxin_free_null_s{seed}"
        runs.append(
            (
                run_id,
                config(
                    run_id=run_id,
                    arm="toxin_free_null",
                    seed=seed,
                    amplitude=null_amplitude,
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
                "arm": cfg["_assay"]["arm"],
                "seed": cfg["_assay"]["seed"],
                "amplitude": cfg["_assay"]["amplitude"],
                "input_relpath": input_path.relative_to(ROOT).as_posix(),
                "input_sha256": digest(input_path),
                "output_relpath": f"generated/results/{index}/output.h5.gz",
            }
        )

    if len(entries) != CONTRACT["design"]["jobs"]:
        raise SystemExit(
            f"REFUSED: generated {len(entries)} runs but the contract declares "
            f"{CONTRACT['design']['jobs']}"
        )

    dump(
        generated / "manifest.json",
        {
            "schema_version": 2,
            "assay_id": CONTRACT["assay_id"],
            "assay_contract_sha256": digest(ROOT / ASSAY_CONTRACT_FILE),
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
