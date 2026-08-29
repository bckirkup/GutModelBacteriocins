# Sealed two-wall Neumann low-screening envelope

The sealed two-wall image series now uses its analytic shell budget for every
positive `kH`, clamped at 512 shells, while retaining the per-shell relative
convergence break. An explicit eight-shell floor is reserved for the
low-screening regime below the cap crossover at `kH = 0.0225`.

The authored accuracy envelope is:

| `kH` | shipped shells | floored | converged at | shipped relative error |
|---:|---:|:---:|---:|---:|
| 0.500000 | 24 | no | 19 | 0 |
| 0.200000 | 58 | no | 45 | 0 |
| 0.158114 | 73 | no | 56 | 0 |
| 0.100000 | 116 | no | 86 | 0 |
| 0.060000 | 192 | no | 137 | 0 |
| 0.050000 | 231 | no | 163 | 0 |
| 0.049900 | 8 | yes | 163 | -5.15e-2 |
| 0.020000 | 8 | yes | 381 | -1.26e-1 |
| 0.005000 | 8 | yes | 1375 | -2.52e-1 |
| 0.001000 | 8 | yes | 6038 | -3.72e-1 |

The dominant measured low-screening error is sealed image-base truncation and
cancellation, not the Robin correction. As `kH` approaches zero, the sealed
problem has no bounded solution: with zero decay and zero flow, mass
accumulates without a sink. That limit needs a policy gate, not more shells.

`qssa.low_screening_policy` is evaluated once during QSSA initialization:

* `warn` (the default) proceeds and reports the measured error scale and
  encountered `kH`;
* `error` rejects the run with `SimulationError`;
* `allow` proceeds silently for deliberate low-screening diagnostics.

Run provenance distinguishes genuine image-series non-convergence
(`neumann_image_series_cap_hits`) from evaluations using the explicit
low-screening floor (`neumann_low_screening_evaluations`). It also records the
pre-clip negative-field count and most-negative value. The final field remains
clamped at zero; the diagnostic is recorded before that clamp.
The negative pre-clip test intentionally uses an unphysical negative source
rate to drive the field negative, so it covers diagnostic counter plumbing
rather than a physically reachable case.

The translated/reflected image construction retains the inherited
wall-normal-flow approximation. This envelope change does not alter that
scientific approximation or the shipped ColE1 `kH = 0.158114` behavior.
