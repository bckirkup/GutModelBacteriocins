# Wall-normal drift and the two-wall Neumann image series

The sealed QSSA bacteriocin kernel sums translated and `z`-reflected images
between the epithelium (`z_lo`) and the lumen (`z_hi`). This document records
what the construction actually solves when the flow has a wall-normal
component, why no per-image reweighting repairs it, the measured error as a
function of the wall-normal Péclet number, and the validity envelope that is
now enforced.

Summary of the result: the defect is real and is exactly first order in the
wall-normal flow, but it is **not** repairable by drift-corrected images — the
exact wall reflection coefficient depends on the transverse wavenumber, and no
real-space image weight can. The series is therefore gated on a measured
envelope rather than "corrected". Shipped defaults sit inside the envelope, and
the shipped numerics are unchanged by this work.

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
∂ψ/∂z + a ψ = 0       at z = z_lo and z = z_hi                      (R)
```

This is the `c_lo = c_hi = +a` sealed reference of `docs/ROBIN_LUMEN_BC.md`,
and it is what `robin_correction_table.cpp` evaluates modally as
`mode_sum(..., robin_boundary = false, ...)`.

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
neighbouring problem with the wrong wall law.

## 3. Why drift-corrected images cannot fix it

The natural repair is to give each image the exponential weight that makes it a
solution of the true drift operator, i.e. work in the gauge variable `ψ` and
reweight the reflections. Solve (I)+(R) exactly instead, in Hankel space:
transforming laterally with wavenumber `q` and writing
`γ = sqrt(q² + κ²)`, the wall reflection coefficient implied by (R) is

```
R(γ) = (γ - a) / (γ + a)
```

and the exact slab solution is the geometric series in `R(γ)²e^{-2γH}` built
from `e^{-γ|z-z_s|}` (this is the closed form used as the second independent
reference in the measurements below).

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

Treatment (a) of the task is therefore rejected on evidence, not on taste.

## 4. Measured error versus `Pe_z`

`Pe_z = U_z H / D_eff`. Reference: the sealed eigenmode expansion (`β_n = nπ/H`
plus the `β = ia` exponential mode, the same construction as
`mode_sum(..., robin_boundary = false, ...)`), converged in mode count and
cross-checked against the independent Hankel closed form of §3 (see the
cross-check paragraph below). Relative field error over nine
source/target/lateral geometries spanning near-wall, mid-slab and long-range
(`ρ` up to 50 µm):

| `Pe_z` | shipped median | shipped max | drift-corrected median | drift-corrected max |
|---:|---:|---:|---:|---:|
| 0 | 3.8e-16 | 1.9e-15 | 3.8e-16 | 1.9e-15 |
| 0.001 | 1.14e-4 | 4.40e-4 | 8.91e-5 | 4.40e-4 |
| 0.010 | 1.14e-3 | 4.39e-3 | 8.90e-4 | 4.39e-3 |
| 0.0463 | 5.26e-3 | 2.01e-2 | 4.09e-3 | 2.01e-2 |
| 0.0926 | 1.05e-2 | 3.95e-2 | 8.10e-3 | 3.95e-2 |
| 0.185 | 2.10e-2 | 7.65e-2 | 1.59e-2 | 7.65e-2 |
| 0.370 | 4.21e-2 | 1.44e-1 | 3.11e-2 | 1.44e-1 |
| 1.50 | 1.65e-1 | 4.20e-1 | 9.31e-2 | 4.20e-1 |
| 6.00 † | 4.69e-1 | 8.98e-1 | 4.12e-1 | 8.98e-1 |
| 25.0 † | 1.31e0 | 5.18e1 | 8.85e-1 | 6.65e1 |

† These rows characterise the image series as mathematics, not the shipped C++
field. Wall-normal flow suppresses the series screening length: with
`kH = (w - |U_z|)H / 2D_eff`, `Pe_z = 1.5` already gives `kH = 0.0187`, below
the `kLowScreeningFloorThreshold = 0.0225` at which the shipped budget floors
to `kLowScreeningShells = 8` (`kH` falls to 4.7e-3 at `Pe_z = 6` and 1.1e-3 at
`Pe_z = 25`). Above `Pe_z ≈ 1.4` the C++ field therefore differs from the
converged series by the low-screening floor as well as by the drift defect, and
the two error terms are not separable there; a direct C++/reference comparison
at `Pe_z = 25` differs from the table by ~0.5 relative for this reason. The
practical consequence is that the two envelopes are *ordered*: the drift gate
at `|Pe_z| = 0.05` fires roughly 28x earlier than the low-screening floor, so a
configuration is warned about drift long before the shell budget degrades. The
rows at and below `Pe_z = 0.37` are unaffected — `kH = 0.065` there — and the
C++ probe reproduces the Python model to near machine precision across the
entire gated range.

The `Pe_z = 0` row is the exactness check: the series reproduces the mode sum
to machine precision, so the measured error is the drift defect and not
truncation. The error is linear in `Pe_z` over the whole low-`Pe_z` range:

```
median relative error ≈ 0.114 · |Pe_z|
max    relative error ≈ 0.439 · |Pe_z|          (|Pe_z| ≲ 0.05)
```

Reversing the flow direction (flow toward the epithelium) is symmetric at small
`|Pe_z|` and worse at large `|Pe_z|` — the wall the advection points into is
where the missing advective flux accumulates.

**Reference cross-check.** Probe by probe, the mode sum and the Hankel form of
§3 agree to ~1e-15 at six of the nine probes. The three that disagree
(1.6e-3 .. 3.4e-3) are exactly those with `z_t = z_s`, where the Hankel
integrand's leading term decays only as `1/γ` and the `J₀` tail is
conditionally convergent. Refining that quadrature drives the disagreement down
as `q_max^{-1/2}`, monotonically toward the mode-sum value (1.3e-2 at
`q_max ρ = 6e2`, 6.6e-4 at `2.5e5`), which identifies it as tail truncation in
the cross-check oracle rather than error in the reference. The mode sum is
separately converged in mode count to 1e-9 by doubling.

**On the audit's residual numbers.** A finite-difference interior residual
cannot be used to compare these constructions at all. A second-order
central-difference residual normalized by `λC`, evaluated at the same point and
spacing for the *exact* free-space kernel (whose true residual is zero) and for
the shipped image sum:

| spacing | exact free-space kernel | shipped image sum |
|---|---:|---:|
| `H/500` | 4.73e1 | 4.72e1 |
| `H/2000` | 4.64e1 | 4.66e1 |
| `H/8000` | 4.63e1 | 4.65e1 |

The exact kernel scores the same value as the image sum to within 0.2% at every
spacing: near the source singularity the metric is entirely discretization
error and cannot distinguish an exact solution from the image sum. I could not
reproduce the reported 2.9e-1 versus 1.4e-3 pair, and whatever it measures it
is not the interior residual of the image sum. The exact residual identity in
§2 makes that comparison unnecessary, and the error table above is a field
comparison against a converged reference, so it does not depend on a difference
stencil.

## 5. Consequence at shipped parameters

`AdvectionField::radial_velocity` uses `profile_alpha = 1.5`, so the
wall-normal flow — and hence the local `Pe_z` a source sees — is depth
dependent:

```
Pe_z(z_s) = Pe_z(z_hi) · (z_s/H)^1.5 ,      Pe_z(z_hi) = 0.0926
```

at shipped `H = 100 µm`, `radial_turnover = 5400 s`, ColE1
`D_eff = 2e-11 m²/s`, `decay_rate = 5e-5 s⁻¹`. Measured relative field error
for near-wall and mid-range targets at each source depth:

| `z_s/H` | local `Pe_z` | median error | max error |
|---:|---:|---:|---:|
| 0.02 | 2.62e-4 | 2.24e-5 | 6.43e-5 |
| 0.05 | 1.04e-3 | 8.45e-5 | 2.47e-4 |
| 0.10 | 2.93e-3 | 2.07e-4 | 6.35e-4 |
| 0.20 | 8.28e-3 | 4.13e-4 | 1.34e-3 |
| 0.30 | 1.52e-2 | 4.85e-4 | 1.60e-3 |
| 0.50 | 3.27e-2 | 1.75e-5 | 2.28e-4 |
| 0.75 | 6.01e-2 | 2.45e-3 | 8.61e-3 |
| 1.00 | 9.26e-2 | 7.13e-3 | 1.77e-2 |

The shipped default founder band is `z/H ∈ [0, 0.5]`
(`Simulation::init_population`), where the measured error is
**6.4e-5 .. 1.6e-3** — consistent with the audit's reported 1.2e-3..6.9e-3
band, at the low end of it, because the depth profile keeps founders in the
weak-flow region. Only sources at the lumen surface reach ~1.8e-2.

The Robin correction does not absorb any of this. `normalized_correction`
returns `robin_modal - sealed_modal` and the caller adds it to the
*image-series* sealed base, so the image error passes through the Robin path
untouched — which is why it dominates that accuracy budget. `toxin.ln_length`
is infinite (Robin disabled) by default, so under shipped defaults the sealed
image series *is* the field.

## 6. The validity envelope, and why there is no drop-in alternative solve

Envelope thresholds read off the measured curve (worst case over the probe
set):

| accuracy target | validated `Pe_z` magnitude |
|---|---:|
| ≤ 1% worst case | 0.023 |
| ≤ 2% worst case | 0.050 |
| ≤ 3% worst case | 0.068 |
| ≤ 5% worst case | 0.114 |

`kDriftEnvelopePeZ = 0.05` (≤ ~2% worst case, ≤ 0.6% median) is the enforced
envelope. `qssa.drift_envelope_policy` mirrors
`qssa.low_screening_policy`: `warn` (default), `error`, `allow`, evaluated once
during QSSA initialization from the same mid-domain flow probe. Run provenance
records `neumann_drift_envelope_evaluations`, the number of kernel evaluations
whose source-local `|Pe_z|` exceeded the envelope, separately from
`neumann_image_series_cap_hits` and
`neumann_low_screening_evaluations`. At shipped defaults the mid-domain probe
gives `Pe_z = 0.033`, inside the envelope, so no warning fires and the counter
stays at zero.

Routing out-of-envelope configurations to an exact alternative solve is *not*
offered, for a measured reason rather than a scheduling one: the only exact
representations available are the modal sum and the Hankel form, and both need
`K₀` per mode per evaluation. `std::cyl_bessel_k` is host-only, so putting
either in the QSSA hot loop would break the shared host/device series that
PR #359 established, and it costs `O(modes)` Bessel evaluations per pair
against the current handful of exponentials. The two routes that could work are
recorded here as follow-up scientific decisions, not implemented:

1. **Drift-consistency table.** Tabulate `robin_modal - sealed_image` on the
   existing Robin correction-table grid instead of
   `robin_modal - sealed_modal`. Then `sealed_image + correction = robin_modal`
   exactly, at zero additional hot-loop cost. This repairs the field only when
   the Robin path is enabled, which is not the default, and it moves shipped
   Robin numbers by up to the §5 error — a scientific decision.
2. **Device-safe `K₀`.** A `__host__ __device__` `K₀` plus a truncated modal
   path would allow a genuine exact route on both backends, at `O(modes)` cost.

## 7. What this change does and does not alter

The image-series arithmetic is untouched: no term, weight, shell budget, or
convergence test changes, so host and device remain bit-identical and every
shipped number is unchanged. What is added is the measured envelope, the policy
gate, the provenance counter, and this derivation.

Device evidence for the CUDA compile comes only from CI: there is no CUDA
toolkit or GPU on the development machine used for this work, and the shared
header change is a `constexpr` and a small `__host__ __device__` predicate.
