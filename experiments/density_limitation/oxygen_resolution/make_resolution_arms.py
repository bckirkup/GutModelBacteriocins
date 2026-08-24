#!/usr/bin/env python3
"""Generate the CPU oxygen resolution-and-timestep measurement arms.

This is a measurement, not a calibration or fitting exercise.  The flora
oxygen sink is disabled in every arm so that the resolution comparison is an
attribution control rather than a competition experiment.  The four completed
arms and the intentionally abandoned one are:

* ``res2_dt60`` -- 2 um cells, 60 s biological step
* ``res4_dt60`` -- 4 um cells, 60 s biological step
* ``res6_dt60`` -- 6 um cells, 60 s biological step
* ``res2_dt10`` -- 2 um cells, 10 s biological step
* ``res2_dt1`` -- 2 um cells, 1 s biological step

Only the horizon, console output cadence, summary cadence, grid spacing, and
biological timestep are changed from the supplied base configuration.
"""

import argparse
import json
from pathlib import Path

BASE_CONFIG = Path("/home/ubuntu/gutibm-campaign/o2compete-oxy/s0_amf1/input.json")
DOMAIN_EXTENTS = (96.0e-6, 96.0e-6, 300.0e-6)
ARMS = {
    "res2_dt60": {"grid_dx": 2.0e-6, "bio_dt": 60.0},
    "res4_dt60": {"grid_dx": 4.0e-6, "bio_dt": 60.0},
    "res6_dt60": {"grid_dx": 6.0e-6, "bio_dt": 60.0},
    "res2_dt10": {"grid_dx": 2.0e-6, "bio_dt": 10.0},
    "res2_dt1": {"grid_dx": 2.0e-6, "bio_dt": 1.0},
}


def _cell_counts(grid_dx: float) -> list[int]:
    counts = [round(extent / grid_dx) for extent in DOMAIN_EXTENTS]
    if any(abs(extent / grid_dx - count) > 1.0e-9
           for extent, count in zip(DOMAIN_EXTENTS, counts)):
        raise ValueError(
            f"grid_dx={grid_dx:g} does not divide every domain extent"
        )
    return counts


def _comments(name: str, grid_dx: float, bio_dt: float,
              counts: list[int]) -> list[str]:
    return [
        "Oxygen resolution-and-timestep measurement; no parameter is fitted.",
        (
            "The flora oxygen sink is off (oxygen.vbf_sink = 0) as the "
            "attribution control."
        ),
        (f"arm = {name}; grid_dx = {grid_dx:.3e} m; bio_dt = {bio_dt:g} s; "
         f"grid = {counts[0]} x {counts[1]} x {counts[2]}."),
        ("Only total_time, output_interval, hdf5.schedule.summary, grid_dx, "
         "and bio_dt differ from the copied base configuration."),
        ("Read funded versus demanded oxygen, maintenance, realized shares, "
         "clips, population, fermentation, oxygen profile, and termination "
         "provenance. Results are not fitted or used to alter parameters."),
    ]


def _arm_config(base: dict, name: str, settings: dict[str, float]) -> dict:
    grid_dx = settings["grid_dx"]
    bio_dt = settings["bio_dt"]
    counts = _cell_counts(grid_dx)
    config = dict(base)
    config["total_time"] = 3600
    config["output_interval"] = 60
    config["hdf5"] = dict(config["hdf5"])
    config["hdf5"]["schedule"] = dict(config["hdf5"]["schedule"])
    config["hdf5"]["schedule"]["summary"] = 1
    config["grid_dx"] = grid_dx
    config["bio_dt"] = bio_dt
    config["_comment"] = _comments(name, grid_dx, bio_dt, counts)
    return config


def main(base_config: Path, output_root: Path) -> None:
    base = json.loads(base_config.read_text())
    output_root.mkdir(parents=True, exist_ok=True)
    for name, settings in ARMS.items():
        config = _arm_config(base, name, settings)
        arm_dir = output_root / name
        arm_dir.mkdir(parents=True, exist_ok=True)
        (arm_dir / "input.json").write_text(
            json.dumps(config, indent=1) + "\n"
        )
        counts = _cell_counts(settings["grid_dx"])
        print(
            f"{name}: grid={counts[0]}x{counts[1]}x{counts[2]} "
            f"grid_dx={settings['grid_dx']:.3e} "
            f"bio_dt={settings['bio_dt']:g}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, default=BASE_CONFIG)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("/home/ubuntu/gutibm-campaign/o2fund-res"),
    )
    arguments = parser.parse_args()
    main(arguments.base_config.expanduser().resolve(),
         arguments.output_root.expanduser().resolve())
