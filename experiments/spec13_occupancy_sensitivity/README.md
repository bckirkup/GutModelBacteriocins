# Spec 13 Layer-3 Mucosal Occupancy Sensitivity (`GAP-05` / Issue `#467`)

This directory stages the compact `2 x 3` factorial sensitivity design evaluating how luminal-to-mucosal reattachment (`layer3.reattach_p0`) interacts with axial luminal carrying profiles (`gradient` vs `uniform`) across seeds `1001–1003` (`18` configurations total).

## Factorial Axes

| Axis | Levels | Parameter Mapping |
| :--- | :--- | :--- |
| Axial colonic profile | `gradient`, `uniform` | `layer3.uniform_profile = false` vs `true` |
| Reattachment probability | `0.1x`, `0.3x`, `1.0x` baseline | `layer3.reattach_p0` in `{0.005, 0.015, 0.05}` |
| Replicate seeds | `1001, 1002, 1003` | `seed` |

## Evaluation Observables (No Arbitrary `< 100%` Binary Gate)

Per the `GAP-05` calibration decision (`docs/CYCLE2_EMPIRICAL_CALIBRATION_REVIEW.md`), `>= 1` cell in a `60 x 60 um` patch is not equivalent to optical microscopy field coverage. Readouts via `analysis/spec13_gradient/run_readout.py` report:
- Strictly post-purge (`t > 48 h`, `173,400–258,600 s`) regional binary occupancy (`N >= 1`), `mean_cells_per_occupied_patch`, and count-thresholded occupancy (`N >= 10`, `N >= 50` cells/patch).
- Last time-series sample (`71h50m`) vs exact `72 h` (`259,200 s`) `prov.json` terminal endpoints.
- Stool `CFU/g` (`1e6–1e8` target band), mucosal flatness (`<= 4.0`), and closed stock ledgers.

## Regenerate Configurations

```bash
python3 experiments/spec13_occupancy_sensitivity/prepare.py
```
