# Wall-normal drift and the two-wall Neumann image series

The sealed QSSA bacteriocin kernel sums translated and `z`-reflected images
between the epithelium (`z_lo`) and the lumen (`z_hi`). This document records
what the construction actually solves when the flow has a wall-normal
component, why no per-image reweighting repairs it, the measured error as a
function of the wall-normal Péclet number, and the validity envelope that is
now enforced.

Summary of the result: the legacy reflected-image series enforces the wrong
sealed wall law when wall-normal flow is present. The physical drift-corrected
path (`qssa.drift_correction=true` or per-source automatic dispatch via
`qssa.drift_envelope_policy=auto` above `qssa.peclet_z_threshold`, default
`0.05`) now evaluates an exact three-branch hybrid Sturm-Liouville / Hankel
sealed kernel directly on both CPU and CUDA GPU (`robin::evaluate_physical_sealed_concentration`),
replacing the earlier 33³ log-rho interpolation table for sealed drift-corrected
sources and eliminating the host-only GPU fallback. On the 300 off-node
geometries at the configured ColE1 basis (`Pe_z = 0.822`), the previous 33³
table's median `1.22e-4`, p90 `3.88e-3`, and max `1.55e-1` error drops to
median `1.65e-14` and max `2.86e-13`. Against both a 500-mode Sturm-Liouville
expansion and a fine-grid axisymmetric finite-volume Robin advection-diffusion
reference across `ColE1` (`R = 50.2`, `Pe_z = 0.822` and `2.325`) and `ColE2`
(`R = 2.04` and `1.2`), relative interior and wall L2 errors are at most
`2.50e-12` (vs Sturm-Liouville) and `1.29e-3` (vs FVM), with zero negative
overshoot and sealed wall flux `|J_z| / (Q / A) <= 1.17e-10`. The shipped
default policy remains `warn` with `qssa.drift_correction=false` to preserve
historical runs bit-for-bit.

## 1. Operator, gauge transform, and the sealed wall law

The QSSA field for one bacteriocin species obeys, with a uniform flow `U`
evaluated at the source (`greens_function.cpp`, `single_kernel`):

```
D ∇²C - U·∇C - λ C + Q δ(x - x_s) = 0
```

`D = D_eff = diff_coeff / retardation`. A genuine sealed wall passes zero
**total** flux, diffusive plus advective:

```
J_z = -D ∂C/∂z + U_z C = 0        at z = z_lo and z = z_hi          (W)
```

Write `a = U_z / (2D)` and gauge out the drift with

```
C = exp(U·(x - x_s) / (2D)) ψ ,   κ² = λ/D + |U|²/(4D²)
```

which turns the interior operator isotropic and the wall law into a Robin
condition with the *same* coefficient at both walls:

```
∇²ψ - κ² ψ = -(Q/D) δ(x - x_s)                                      (I)
∂ψ/∂z - a ψ = 0       at z = z_lo and z = z_hi                      (R)
```

This is the physical `eigen_a = -a` sealed reference evaluated by
`robin_correction_table.cpp` as `mode_sum(..., robin_boundary = false, ...)`.
Its zero mode is `exp(+a zeta)`, so the gauge-transformed field carries the
physical equilibrium profile. The legacy image-consistent reference retains
`eigen_a = +a` only for the default-off Robin-table subtrahend.

For `a = 0`, (R) is `∂ψ/∂z = 0`: plain mirror images are exact, and the shipped
series is exact. Everything below is about `a ≠ 0`.

## 2. What the shipped series actually solves

The shipped reflected image is the free-space kernel placed at the mirrored
source with the wall-normal flow component reversed
(`concentration_sealed`: `image_flow[2] = -image_flow[2]`). Two consequences,
both exact:

**(a) Interior residual, first order in `U_z`.** A reflected image `f` solves
the operator built from `U' = (U_x, U_y, -U_z)`, so applying the *physical*
operator `L = D∇² - U·∇ - λ` to it leaves

```
L f = (U' - U)·∇f = -2 U_z ∂f/∂z
```

The image sum therefore has interior residual `-2 U_z Σ_reflected ∂f_i/∂z`: not
a solution of the physical PDE, not a solution of the mirrored one, and exactly
proportional to the wall-normal velocity. It vanishes identically for
wall-parallel flow.

**(b) The wall condition enforced is the wrong one.** Reversing `U_z` makes the
reflected image the exact mirror of the direct kernel about `z_lo`, so the pair
satisfies `∂C/∂z = 0` at the wall — zero *diffusive* flux. The physical sealed
wall is (W). The shipped construction therefore leaks advective flux `U_z C` at
both walls; a sealed slab built this way does not conserve mass when `U_z ≠ 0`.

Statement (b) is the more useful description of the defect: the series is not
"an approximate solution of the right problem", it is the exact solution of a
neighbouring problem with the wrong wall law. The previous in-tree sealed modal
field used the same non-physical zero-diffusive-flux wall law as the image
series: its wall-flux residual was approximately `+0.451`, whereas the physical
`eigen_a = -a` modal field gives approximately `2e-10`.

## 3. Why per-image reweighting cannot fix it, and the table correction

The natural repair is to give each image the exponential weight that makes it a
solution of the true drift operator, i.e. work in the gauge variable `ψ` and
reweight the reflections. Solve (I)+(R) exactly instead, in Hankel space:
transforming laterally with wavenumber `q` and writing
`γ = sqrt(q² + κ²)`, the wall reflection coefficient implied by (R) is

```
R(γ) = (γ - a) / (γ + a)
```

and the exact slab solution is the geometric series in `R(γ)²e^{-2γH}` built
from `e^{-γ|z-z_s|}`. This is the closed form used as a reference in the
measurements below, but not an independent check of the legacy image family:
it uses the same reflection coefficient implied by that family.

`R` depends on `γ`, hence on the transverse wavenumber `q`. A real-space image
carrying a constant multiplicative weight — which is what
`exp(U_z (z - z_image) / (2D))` and every other per-image exponential factor
is, since it does not depend on the lateral separation — has a `q`-independent
reflection coefficient. **No constant-weight image family can satisfy (R) at
both walls for `a ≠ 0`.** The construction is exact only in the limits `a → 0`
and `γ ≫ a` (short range, where `R → 1` and plain mirror images are recovered),
which is precisely why the measured error below is smallest for near-wall,
short-range probes and largest across the slab.

This is not a statement about one candidate weight. Scanning *all* constant
reflection weights `w` (the shipped series is `w = 1`; every per-image
exponential correction is some member of this family) and taking the best
achievable worst-case error over the probe set gives:

| `Pe_z` | shipped (`w = 1`) | best `w` | best achievable | gain |
|---:|---:|---:|---:|---:|
| 0.0463 | 2.01e-2 | 1.0135 | 1.35e-2 | 1.48x |
| 0.0926 | 3.95e-2 | 1.0264 | 2.69e-2 | 1.47x |
| 0.1850 | 7.65e-2 | 1.0508 | 5.32e-2 | 1.44x |
| 0.3700 | 1.44e-1 | 1.0939 | 1.04e-1 | 1.38x |
| 1.5000 | 4.20e-1 | 1.2172 | 3.58e-1 | 1.17x |

The entire family buys at most a factor ~1.5 and does not change the scaling:
the error stays first order in `Pe_z`. The specific
`exp(U_z(z - z_image)/(2D))` weighting was also measured directly and improves
the *median* error by ~23%.

Its *worst case* is unchanged to four digits, and that deserves an explanation
rather than a footnote, because it is a property of the probe set and not of
the correction. The worst case over the nine probes is attained at
`(z_s, z_t) = (0.05H, 0.95H)`, and at every probe with `z_s + z_t = H` exactly
the reversed-`U_z` and non-reversed reflected families coincide term by term,
so any change confined to the reflected images is invisible there. Measured at
an asymmetric probe (`z_s = 0.40H`, `z_t = 0.98H`, `ρ = 5 µm`) the two
constructions separate cleanly — shipped 4.27e-3 versus 1.03e-2 without the
reversal at `Pe_z = 0.0926`, a factor 2.4 — which is the sensitivity the
regression test asserts, precisely so that a change to the reflected family
cannot pass unnoticed behind a symmetric worst case.

Per-image reweighting is therefore rejected on evidence, not on taste. The
initial PR #391 implementation tabulated
`T_drift = sealed_modal(physical) - image_series(runtime)` on a 33³ log-rho
grid and fell back to the host on GPU. That table still exhibited up to
`1.55e-1` off-node relative error for three reasons: (1) near the walls,
`T_drift` retained a singular `(1 - exp(-2a ζ)) / r_lo` reflection mismatch
that trilinear interpolation on 33 nodes could not resolve; (2) linear
interpolation of the exponentially screened field `exp(-κ ρ)` across 33
log-spaced nodes (where `Δρ ~ 35 µm` near the cutoff against `1/κ ~ 3.9 µm`
for retarded ColE1 under distal flow) incurred severe interpolation error; and
(3) the near-wall direct path previously used a radial-only gauge factor rather
than the 3D angle-aware gauge `exp((U_x Δx + U_y Δy + U_z Δz) / (2D))`.

Issue #478 replaces sealed drift-corrected table lookup with a shared
`__host__ __device__` three-branch hybrid Sturm-Liouville / Hankel kernel
(`robin::evaluate_physical_sealed_concentration` in
`src/diffusion/robin_correction_table.h`) based on the exact decomposition
`ψ_phys(ρ, z_s, z_t) = ψ_N(ρ, z_s, z_t; κ) + F(ρ, z_s + z_t)`, where `ψ_N` is
the pure-Neumann image series (satisfying `∂ψ/∂z = 0`) and `F` depends only on
`σ = z_s + z_t` and `ρ`:

- **Strong screening (`(κ - |a|) H >= 2.0`, e.g., retarded ColE1 under distal
  flow):** pure-Neumann image series plus non-oscillatory positive real-space
  drift integrals `I_+` and `I_-` evaluated with 24-point Gauss-Legendre
  quadrature with the 3D exponential gauge folded inside the integrand to avoid
  overflow and catastrophic cancellation.
- **Weak screening, mid/far lateral range (`(κ - |a|) H < 2.0`, `ρ >= 0.15 H`):**
  exact closed-form physical Sturm-Liouville eigenmode sum (`β_n = nπ/H`,
  scaled `K_0e` with monotonic early exit).
- **Weak screening, near-axis (`(κ - |a|) H < 2.0`, `ρ < 0.15 H`, including
  `ρ = 0`):** pure-Neumann image series plus exact `m = 0` `I_±` wall integrals
  and a rapidly decaying `exp(-2γH)` Hankel tail with pole-cancelling
  substitution `u = δ (exp(s) - 1)`.

For finite `lumen_transfer_length` Robin boundaries with drift correction,
`normalized_robin_radial` restores the full 3D angle-aware gauge on the direct
near-wall path while the smooth Robin-minus-physical-sealed correction table
retains `table->physical_values`.

## 4. Measured error versus `Pe_z`

The physical reference uses `eigen_a = -a`; the legacy image-consistent modal
field uses `eigen_a = +a`. The corrected table uses the physical reference and
subtracts the same image series that the runtime evaluates, including its
relative tolerance, shell budget, explicit-shell setting, and legacy-reflection
mode. The settled measurements against the physical reference are:

| `Pe_z` | image vs physical, median | max | image vs legacy `+a`, median | max |
|---:|---:|---:|---:|---:|
| 0 | `1.9e-10` | `2.3e-10` | `1.9e-10` | `2.3e-10` |
| 0.001 | `1.67e-4` | `4.40e-4` | `1.14e-4` | `4.40e-4` |
| 0.010 | `1.67e-3` | `4.39e-3` | `1.14e-3` | `4.39e-3` |
| 0.0463 | `7.70e-3` | `2.01e-2` | `5.26e-3` | `2.01e-2` |
| 0.0926 | `1.54e-2` | `3.95e-2` | `1.05e-2` | `3.95e-2` |
| 0.185 | `3.06e-2` | `7.65e-2` | `2.10e-2` | `7.65e-2` |
| 0.370 | `6.07e-2` | `1.44e-1` | `4.20e-2` | `1.44e-1` |
| 0.822 | `1.32e-1` | `2.78e-1` | `9.27e-2` | `2.78e-1` |
| 1.50 | `2.53e-1` | `5.17e-1` | `1.80e-1` | `5.17e-1` |

The corrected fitted envelope for `|Pe_z| <~ 0.05` is therefore:

```
median relative error ~ 0.167 |Pe_z|
max    relative error ~ 0.439 |Pe_z|
```
The earlier `0.114*|Pe_z|` and `0.439*|Pe_z|` coefficients were calibrated
against the wrong legacy `+a` modal reference. The physically referenced
coefficients above are the settled nine-probe values. The earlier Hankel
cross-check was not independent evidence: it encoded the same reflection
coefficient as the legacy image family.

The legacy image arithmetic remains unchanged when the correction is disabled.
When enabled (`qssa.drift_correction=true` or per-source
`qssa.drift_envelope_policy=auto`), the sealed path evaluates the exact
three-branch hybrid Sturm-Liouville / Hankel kernel directly on both CPU and
CUDA GPU without table interpolation. On the 300 off-node geometries at the
configured ColE1 basis (`Pe_z = 0.822`), the residual against the converged
physical Sturm-Liouville reference drops from the PR #391 33³ log-rho table's
median `1.22e-4` / p90 `3.88e-3` / max `1.55e-1` to median `1.65e-14` / max
`2.86e-13`.

Measured before/after accuracy across `ColE1` (`R = 50.2`) and `ColE2`
(`R = 2.04` and `1.2`) over `z_s in {0.1H, 0.5H, 0.9H}` (`tests/test_wall_normal_peclet.cpp`,
zero negative overshoot in all cases):

| Benchmark case | `Pe_z` | Legacy image wall L2 (vs 500-mode SL) | Hybrid interior L2 (vs 500-mode SL) | Hybrid wall L2 (vs 500-mode SL) | Hybrid interior L2 (vs FVM) | Hybrid wall L2 (vs FVM) | Max sealed wall flux `\|J_z\|/(Q/A)` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ColE1` mid-depth (`R = 50.2`) | 0.822 | `2.37e-1` | `1.27e-13` | `1.20e-13` | `1.07e-3` | `1.28e-3` | `4.06e-11` |
| `ColE1` near-lumen (`R = 50.2`) | 2.324 | `5.02e-1` | `1.78e-14` | `3.01e-13` | `1.13e-3` | `1.29e-3` | `5.02e-11` |
| `ColE2` near-lumen (`R = 2.04`) | 0.108 | `4.39e-2` | `1.83e-12` | `1.49e-12` | `9.08e-4` | `9.69e-4` | `1.03e-10` |
| `ColE2` near-lumen (`R = 1.2`) | 0.064 | `2.64e-2` | `2.18e-12` | `2.50e-12` | `5.20e-4` | `5.12e-4` | `1.17e-10` |
| 300 off-node geometries (`ColE1`) | 0.822 | `1.55e-1` max (PR #391 33³ table) | — | `2.86e-13` max (`1.65e-14` median) | — | — | — |

Under combined 3D distal plus wall-normal flow (`U_x != 0`, `U_z != 0`,
`Pe_z = 0.822` on `ColE1`), the hybrid kernel satisfies the physical sealed
wall law with max relative wall residual `4.59e-8` across `z_lo` and `z_hi`.

The Hankel cross-check is useful for the physical Robin mathematics, but it is
not independent evidence against the legacy image family: it is built from
`R(γ) = (γ-a)/(γ+a)`, the reflection coefficient implied by that family. The
agreement to approximately `1e-15` is consequently expected.

## 5. Consequence at shipped parameters

`AdvectionField::radial_velocity` uses `profile_alpha = 1.5`, so the
wall-normal flow — and hence the local `Pe_z` a source sees — is depth
dependent. For `H = 1e-4 m` and `radial_turnover = 5400 s`,
`U_z(0.5H) = 6.547e-9 m/s` and `U_z(H) = 1.8519e-8 m/s`.

```
U_z(z) = (H / radial_turnover) · (z/H)^1.5
```

The initialization gate's configured ChemicalSpec basis uses
`diff_coeff = 4e-11 m²/s`, `retardation = 10`, hence
`D_eff = 4e-12 m²/s` and midpoint-probe `Pe_z = 0.164`, outside the enforced
`0.05` envelope. Runtime Green's-function sources instead use each plasmid BI
locus's `diff_coeff` and pI-derived retardation. The configured ColE1 basis has
`Pe_z = 0.822` at mid-depth and `Pe_z = 2.325` at the lumen surface. Thus the
default warning fires, and the runtime basis can be substantially more
restrictive than the ChemicalSpec basis.

The source-local `Pe_z` values remain depth dependent. The corrected physical
reference uses the settled envelope in §4; the previous species table's median
and worst-error columns were based on the legacy `+a` reference and are not
reused as physical claims here.

The strongly pI-retarded ColE1 and colicin M therefore have measured errors of
**5–15% median and 20–40% worst case**, while the audit's sub-percent
characterization applies only to weakly retarded toxins such as ColB, colicin
E2, and microcin V.

When correction is disabled, the Robin correction retains the legacy composition:
`normalized_correction` returns `robin_modal - sealed_modal` and the caller adds
it to the image-series sealed base, preserving shipped values. With correction
enabled, the sealed base uses `evaluate_physical_sealed_concentration` directly
and the smooth Robin table stores `robin_modal - physical_sealed_modal`; near-wall
Robin evaluations use `normalized_robin_radial` with the full 3D angle-aware
gauge and do not add the image base twice. `toxin.ln_length` is infinite (Robin
disabled) by default, so under shipped defaults the sealed image series remains
the field.

## 6. The validity envelope and runtime policy

Envelope thresholds read off the measured curve (worst case over the probe
set):

| accuracy target | validated `Pe_z` magnitude |
|---|---:|
| ≤ 1% worst case | 0.023 |
| ≤ 2% worst case | 0.050 |
| ≤ 3% worst case | 0.068 |
| ≤ 5% worst case | 0.114 |

`qssa.peclet_z_threshold` (default `kDriftEnvelopePeZ = 0.05`, ≤ ~2% worst case
for the legacy image series) configures the wall-normal Péclet threshold.
`qssa.drift_envelope_policy` accepts four values:

- `warn` (default): logs a warning during QSSA initialization whenever the
  configured-species or configured-plasmid midpoint `|Pe_z|` exceeds
  `qssa.peclet_z_threshold`, and retains the legacy image series unless
  `qssa.drift_correction=true`.
- `auto`: automatically promotes each source whose local `|U_z| > 0` and
  `|Pe_z| > qssa.peclet_z_threshold` to the exact hybrid physical kernel while
  keeping low-`Pe_z` sources on the fast image series; logs a concise `NOTE`
  during initialization rather than a warning.
- `error`: fails closed (throwing `SimulationError` naming the offending
  species/plasmid and `Pe_z`) at initialization and whenever a runtime source
  exceeds `qssa.peclet_z_threshold` without drift correction enabled.
- `allow`: suppresses the warning and retains the legacy image series for
  deliberate diagnostic comparisons.

Run provenance under `/run_provenance/` records:

- `neumann_drift_correction`: whether global `qssa.drift_correction` was enabled.
- `qssa_drift_envelope_policy`: configured policy (`warn`, `auto`, `error`, `allow`).
- `qssa_peclet_z_threshold`: configured `|Pe_z|` threshold.
- `neumann_drift_envelope_evaluations`: kernel evaluations whose source-local
  `|Pe_z|` exceeded the threshold (separate from `neumann_image_series_cap_hits`
  and `neumann_low_screening_evaluations`).
- `drift_modal_evaluations`: kernel evaluations executed by the physical
  Sturm-Liouville / Hankel drift kernel (summed across MPI ranks and folded
  from device atomics on GPU).
- `peclet_auto_dispatched_sources`: source solves automatically promoted to the
  physical drift kernel under `qssa.drift_envelope_policy=auto` (counted once
  per owned source before MPI broadcast).

## 7. What this change does and does not alter

With `qssa.drift_correction` unset or false and `qssa.drift_envelope_policy=warn`
(shipped defaults), the sealed and bounded paths retain the legacy
image-consistent table and interpolation axis, preserving shipped results
bit-for-bit within a build. The default-off regression is that same-binary
contract, not a hexadecimal pin: the shipped parameter object must match an
explicit disabled correction, and the legacy median error at the ColE1 mid-depth
basis must stay at least 0.05 (measured 0.151 on the five probes; the §4
nine-probe median is 0.132) rather than collapse onto the physical modal field.
Literals captured with GCC 13 matched Ubuntu 24.04 CI and diverged when
`ubuntu-latest` moved to GCC 15.2, without any change to the series.

With `qssa.drift_correction=true` or `qssa.drift_envelope_policy=auto`, both
CPU and CUDA GPU execute the same shared `__host__ __device__` hybrid physical
kernel (`robin::evaluate_physical_sealed_concentration`) directly in
`src/diffusion/robin_correction_table.h` and `src/gpu/gpu_common.cuh`. Sealed
drift-corrected and auto-dispatched sources no longer fall back to the host on
GPU (`robin_host_fallback_sources == 0`); only near-wall finite-`lumen_transfer_length`
Robin sources retain the existing host fallback. `tests/test_greens_function_gpu.cpp`
asserts `<= 1e-6` relative L2 CPU/GPU parity on both high-`Pe_z` `ColE1`
drift-corrected fields and mixed `ColE1`/`ColE2` auto-dispatched fields, and
verifies zero host fallback on the T4 `gpu-device-tests` CI job.

