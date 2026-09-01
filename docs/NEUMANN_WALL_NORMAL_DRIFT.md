# Wall-normal drift and the two-wall Neumann image series

The sealed QSSA bacteriocin kernel sums translated and `z`-reflected images
between the epithelium (`z_lo`) and the lumen (`z_hi`). This document records
what the construction actually solves when the flow has a wall-normal
component, why no per-image reweighting repairs it, the measured error as a
function of the wall-normal Péclet number, and the validity envelope that is
now enforced.

Summary of the result: the legacy reflected-image series enforces the wrong
sealed wall law when wall-normal flow is present. The opt-in
`qssa.drift_correction` path now uses the physical modal sealed field, a 33³
log-rho correction table against the exact runtime image series, and direct
modal evaluation for the near-wall singular geometry. The correction is off by
default, preserving shipped values bit-for-bit; drift-corrected GPU diffusion
uses the host fallback rather than silently diverging. The measured
interpolation residual is median `1.22e-4`, p90 `3.88e-3`, and max `1.55e-1`
over 300 off-node geometries.

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
implemented correction does not reweight images: it tabulates
`T_drift = sealed_modal(physical) - image_series(runtime)` on the log-rho grid,
then adds that difference to the runtime image base. This preserves the exact
series budget and makes the corrected composition converge to the physical
modal sealed field without a hot-loop Bessel evaluation.

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
When enabled, the correction table uses 33 nodes in each of `(z_source,
z_target,rho)`, with `rho` geometrically spaced from `0.25e-6` to the cutoff
and interpolation performed on `log(rho / rho_min)`. The same log-rho axis is
used for the Robin correction values. At the singular near-wall geometries,
`requires_direct_evaluation()` bypasses interpolation and evaluates the
physical modal field directly. The measured off-node correction residual is
median `1.22e-4`, p90 `3.88e-3`, and max `1.55e-1` over 300 geometries at the
configured ColE1 basis.

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
enabled, the table stores `robin_modal - physical_sealed_modal` and the sealed
base adds the independently tabulated `physical_sealed_modal - image_series`;
the composition therefore approaches the physical Robin modal field. Near-wall
cases use the direct physical modal field and do not add the image base twice.
`toxin.ln_length` is infinite (Robin disabled) by default, so under shipped
defaults the sealed image series remains the field.

## 6. The validity envelope and runtime policy

Envelope thresholds read off the measured curve (worst case over the probe
set):

| accuracy target | validated `Pe_z` magnitude |
|---|---:|
| ≤ 1% worst case | 0.023 |
| ≤ 2% worst case | 0.050 |
| ≤ 3% worst case | 0.068 |
| ≤ 5% worst case | 0.114 |

`kDriftEnvelopePeZ = 0.05` (≤ ~2% worst case, with median approximately
`0.167*|Pe_z|` in the settled low-drift fit) remains the legacy-image
validity envelope. `qssa.drift_envelope_policy` mirrors
`qssa.low_screening_policy`: `warn` (default), `error`, `allow`, evaluated once
during QSSA initialization from the same mid-domain flow probe. Run provenance
records `neumann_drift_envelope_evaluations`, the number of kernel evaluations
whose source-local `|Pe_z|` exceeded the envelope, separately from
`neumann_image_series_cap_hits` and
`neumann_low_screening_evaluations`. At shipped defaults the ChemicalSpec midpoint basis gives `Pe_z = 0.164` and
the configured-plasmid ColE1 basis gives `Pe_z = 0.822`, both outside the
legacy envelope. With `qssa.drift_correction=true`, this Pe_z warning is
suppressed and replaced by a note describing the measured interpolation
residual envelope.

The correction is deliberately opt-in through `qssa.drift_correction` and does
not change any shipped default. Its provenance dataset is
`/run_provenance/neumann_drift_correction`; the existing
`neumann_drift_envelope_evaluations` counter remains available for legacy image
series evaluations. When correction is enabled, the initialization message
states the measured interpolation envelope rather than treating `|Pe_z| > 0.05`
alone as a failure. Device table upload and device-side corrected evaluation are
follow-up work; drift-corrected GPU diffusion is routed to the host fallback.

## 7. What this change does and does not alter

With `qssa.drift_correction` unset or false, the sealed and bounded paths retain
the legacy image-consistent table and interpolation axis, preserving shipped
results bit-for-bit. The opt-in host correction adds the physical sealed modal
reference and the runtime-budget-matched drift table; it records
`/run_provenance/neumann_drift_correction` and preserves the existing
`neumann_drift_envelope_evaluations` counter. Near-wall direct evaluation avoids
adding the image base twice. Drift-corrected GPU diffusion is declined to the
host fallback until device table upload and T4 evidence are available.

Device evidence for CUDA remains a CI concern: there is no CUDA toolkit or GPU
on the development machine used for this work.
