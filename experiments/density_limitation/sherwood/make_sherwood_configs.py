"""Sherwood-vs-none arms: does capping uptake at the diffusive delivery to each
cell create a carrying capacity below the dysbiosis guard?

Lead-authored. Built from the bracket configs so the only differences from the
already-analysed bracket are `uptake_limit` and the arm set. Three flux values
spanning the measured zero-crossing of net growth (0.14x negative, 0.18x/0.22x
positive over the first 17 h), each run with `uptake_limit=none` (the bracket's
behaviour, as the paired control) and `uptake_limit=sherwood`.

kd_b12_btuB is held at the 1e-6 default: the bracket showed the corrinoid
regime does not measurably move guard-crossing time, so it is not a confound
here, and holding it fixed keeps these six arms comparable to the bracket's
kd1e6 column run-for-run.

The run root contains bracket/bracket_f*_kd1e6_s1001.json as source configs;
generated configs are written under sherwood/.
"""

import argparse
import json
import os
import pathlib

MULTIPLIERS = [0.14, 0.18, 0.22]
LIMITS = ["none", "sherwood"]


def main(run_root: pathlib.Path) -> None:
    source_dir = run_root / "bracket"
    out = run_root / "sherwood"
    out.mkdir(parents=True, exist_ok=True)
    for mult in MULTIPLIERS:
        tag = f"f{round(mult * 100):03d}"
        base = json.loads(
            (source_dir / f"bracket_{tag}_kd1e6_s1001.json").read_text()
        )
        for limit in LIMITS:
            cfg = dict(base)
            cfg["uptake_limit"] = limit
            cfg["_comment"] = [
                "Sherwood uptake-limit arm.",
                (f"uptake_limit = {limit}; carbon.epithelial_flux = "
                 f"{mult:.2f}x J_dir ({cfg['carbon.epithelial_flux']:.5e} "
                 "mol/m^2/s); kd_b12_btuB at the 1e-6 default."),
                "Question: the delivery bracket found no flux that holds a",
                "population for 168 h - division rate is linear in delivery",
                "(0.72x mult per agent per hour) while losses are flat at",
                "0.11/agent/h, so net growth crosses zero once with nothing",
                "bending it back. Agents book 4-15% of epithelial carbon and",
                "funded/demanded was exactly 1.0000, i.e. uptake_limit=none",
                "has no density-dependent brake at all. Sherwood caps uptake",
                "at 4*pi*D_eff*r*C_local, which tightens as cells crowd.",
                "Outcome read per arm: does the guard fire, at what time and",
                "density, what fraction of agents are uptake-bound, and how",
                "far funded/demanded falls below 1.",
                "The paired uptake_limit=none arm is the control and also the",
                "sherwood-vs-none divergence measurement.",
            ]
            name = f"sherwood_{tag}_{limit}_s1001.json"
            (out / name).write_text(json.dumps(cfg, indent=1) + "\n")
            print(f"{name}: flux={cfg['carbon.epithelial_flux']:.5e} "
                  f"uptake_limit={limit}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        default=os.environ.get("GUTIBM_CAMPAIGN_ROOT", "."),
        help="campaign root containing bracket/ and sherwood/ (default: "
        "GUTIBM_CAMPAIGN_ROOT or current directory)",
    )
    main(pathlib.Path(parser.parse_args().run_root).expanduser().resolve())
