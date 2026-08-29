# Robin lumen boundary for the QSSA toxin Green's function — design spec

Status: implemented. Numbers below are measured, not assumed; independent
reference calculations are described where they are used.

## 1. What is being replaced and why

The corrected two-wall series from #359 enforces no flux at **both** walls, i.e.
a sealed slab. `AdvectionField::radial_velocity` is nonzero at `z_hi`, so the
lumen wall is not physically sealed, and the choice matters far more than the
bug #359 fixed: measured at shipped parameters (`H=100 µm`, `D_eff=2e-11 m²/s`,
`kH=0.158`), a transparent lumen gives 0.35–0.67× and a perfect sink 0.01–0.52×
the sealed toxin concentration, against a 9–19% bias
from the image-series defect itself.

Advection cannot carry the escape on its own: `Pe = U_z·H/D_eff ≈ 0.086`, so
diffusive escape outruns advective removal by roughly 12:1. The screening length
`1/k ≈ 633 µm` is 6× the slab, so the field shape is set by the lumen boundary
condition, not by decay.

Decision (Benjamin, this session): Robin mass transfer at `z_hi`, no flux at
`z_lo`. Sealed and perfect-sink are its two limits, so the model spans the whole
bracket with one parameter instead of picking an extreme.

## 2. Governing problem

Screened advection–diffusion for one toxin species, QSSA, uniform flow `U`
evaluated at the source (unchanged from current code):

```
-D ∇²C + U·∇C + λC = Q δ(r - r_s)
z = z_lo:  -D ∂C/∂z + U_z C = 0        (epithelium, zero total flux)
z = z_hi:  -D ∂C/∂z = k_c C            (lumen mass transfer)
lateral:   periodic (unchanged)
```

`k_c` [m/s] is the lumen-side mass-transfer coefficient. `Bi = k_c·H/D_eff`.

## 3. Drift transform and why both walls become Robin

Substituting `C = exp(U·(x - x_s)/(2D)) ψ` removes the drift from the interior:

```
-D ∇²ψ + (λ + |U|²/(4D)) ψ = Q δ(r - r_s)
```

The transform is *not* boundary-neutral. Writing `a = U_z/(2D)` and
`c_hi = k_c/D + a`, the two conditions become

```
z_lo:  ψ' - a ψ = 0
z_hi:  ψ' + (k_c/D + a) ψ = 0
```

so the transformed problem is Robin at **both** walls even though the physical
condition at `z_lo` is pure no-flux. Ignoring this would be a `Pe/2 ≈ 0.043`
error in the mode structure at shipped parameters — small, but it costs nothing
to carry, because the eigenvalue solve happens once at init.

Here `c_lo = -a` and `c_hi = k_c/D + a`. The drift sign differs at the two
walls because the epithelial law is zero total flux while the lumen law is
diffusive film transfer on top of free advective outflow.

## 4. Eigenmodes

```
φ_n(z) = cos(β_n ζ) - (a/β_n) sin(β_n ζ),      ζ = z - z_lo
tan(β_n H) = (c_hi - c_lo) β_n / (β_n² + c_lo c_hi),
  c_lo = -a
N_n = ∫₀^H φ_n(ζ)² dζ
κ_n = sqrt( (λ + |U|²/(4D)) / D + β_n² )

C(r) = exp(U·(x - x_s)/(2D)) · (Q / 2πD) · Σ_n φ_n(ζ_s) φ_n(ζ_t) K₀(κ_n ρ) / N_n
```

`ρ` is the lateral (xy) separation. The roots are found by bracket-and-bisect
on each branch, with diagnostics including `c_lo`, `c_hi`, `Bi`, and `aH`.
For these mixed signs, no imaginary root is present and branch zero is
bracketed directly. The sealed reference is handled as a separate mode set;
its exponential mode is described below.

For the sealed reference, `c_lo = c_hi = +a`. When `a != 0`, beta zero is
replaced by the imaginary mode `beta = i*a`, with
`phi = exp(-a*zeta)` and
`N = (1 - exp(-2*a*H))/(2*a)` (using `N=H` when `a=0`). Its radial decay is
`kappa² = lambda/D + (U_x² + U_y²)/(4D²)`. Only the exact `a=0` limit uses
the beta-zero mode.

## 5. Why this is NOT a runtime mode sum

Measured mode counts for 1e-6 convergence at `Bi = 4` (same script):

| ρ | modes |
|---|---|
| 1 µm | 393 |
| 2 µm | 197 |
| 5 µm | 79 |
| 20 µm | 21 |
| 50 µm | 9 |

`K₀` converges slowly near the source, and near-source is where agents evaluate.
This retracts the "10–20 modes, likely cheaper than 73 image shells" estimate I
gave earlier — it is true only in the far field. Summing the *correction* instead
of the field only halves the count (199 at ρ=1 µm), so that is not a way out
either.

## 6. Implementation: sealed images + tabulated correction

Keep #359's image series for the sealed part and add a per-species table of

```
ΔĈ(z_s, z_t, ρ) := (4πD/Q) · [ C_robin_modes - C_sealed_modes ]
C_robin ≈ C_sealed_images + (Q/4πD) · ΔĈ_interp
```

In the interior, `ΔC` is regular at `ρ → 0` because the free-space singularity
is identical in both fields and cancels. When both source and target are within
a cell radius of a wall, near-axis cancellation is not bounded on the table's
z scale; the evaluator therefore routes those queries to the independent
direct mode sum. Elsewhere the correction tabulates on a uniform grid
**including** `ρ = 0`. Runtime cost is one trilinear interpolation, identical
on host and device, with no Bessel evaluation in the hot loop.

- Table axes: `z_s, z_t ∈ [z_lo, z_hi]`, `ρ ∈ [0, cutoff]`, uniform, `n = 33`
  per axis (0.27 MB per species in doubles).
- Built once at init from the mode sum with a high mode count (target 2000+;
  `ΔĈ` is converged by 600 at ρ=0.6 µm).
- One eigenvalue set per `z_s` row: `a` and `c_hi` depend on the source-height
  mean profile `U(z_s)`. The correction table deliberately neglects peristaltic
  modulation so its key remains stationary; the sealed image base and runtime
  drift restoration retain the full time-varying flow.
- Scale by `Q` at runtime; the tabulated quantity is `Q`-independent.

Measured interpolation error against the direct mode sum is at most `1.8e-3`
for shipped-flow near-lumen samples. The accuracy budget is, in descending
order of importance: (1) inherited image-base inconsistency under wall-normal
flow, `1.2e-3..6.9e-3` at shipped `U_z`, (2) rho/z interpolation, `<=1.8e-3`,
and (3) sealed-reference mode-set error, `<=1.4e-3`. Defining the correction
against the image series was rejected: its near-axis difference contains a
singular component and the uniform rho table has `14%-528%` error.

## 7. Parameter: expose the length, not the coefficient

Config key `toxin.lumen_transfer_length` = `δ` [m], with

```
k_c = D_eff / δ         (default effective basis)
k_c = D_free / δ        (optional free basis)
Bi  = k_c·H/D_eff
```

`δ` is the measurable quantity (a lumen-side unstirred/diffusive boundary layer)
and making it the input keeps the two species consistent automatically.

Default `δ = 100 µm` → `Bi = 1` for the effective basis; the optional free
basis gives `Bi = R`. The effective basis is the default because it avoids
making lumen absorption scale with poorly constrained retardation. Basis:
human jejunal unstirred-layer thickness
is 35–40 µm (Levitt et al., *J Clin Invest* 1990, range 23–65; Strocchi &
Levitt, *Am J Physiol* 1992, ≈35 µm), which replaced older ~600 µm osmotic
estimates. Those are small-bowel measurements with glucose probes; the colon is
less stirred and carries an outer mucus layer, so the default is deliberately
2–3× the jejunal value and is documented as **colon-adjusted, not measured**.

`Bi = 0.3 / 2 / 5` (δ ≈ 600 / 100 / 40 µm) puts near-lumen toxin at ≈0.75× /
0.45× / 0.2× the sealed value, so `δ` must appear in the sensitivity sweep for
any lethal-core/halo or comet-tail claim regardless of the default.

At low screening, `|Delta|/C_robin` rises to `2.4` at `kH=5e-3`, a mild
cancellation of approximately `0.4` digits. This is a documented soft spot,
not a gate; the direct mode path remains reachable for near-wall/near-axis
queries.

The direct fallback uses
`rho < cell_radius && source_wall_distance < cell_radius &&
target_wall_distance < cell_radius`. Requiring both points to be near the same
wall keeps this expensive safety-net path bounded while retaining the
near-wall/near-axis protection. Direct-mode evaluations are counted in
`run_provenance/robin_direct_mode_evaluations`, alongside the table build and
eviction counters.

Limits to keep reachable and tested: `δ → ∞` (`k_c = 0`) must return the #359
sealed result bit-for-bit via the same code path, and large `Bi` must drive
`C(z_hi) → 0`.

## 8. Required tests (independent oracles, not self-oracles)

1. **Flux residual at both walls**, numerically differentiated: `-D ∂C/∂z + U_z C = 0`
   at `z_lo` and `-D ∂C/∂z - k_c C = 0` at `z_hi`, over several source
   positions and lateral offsets, normalised by `D·C/H`.
2. **Table vs direct mode sum**: max relative error over shipped-flow near-lumen
   points below 5e-3 of the local total field (measured <=1.8e-3 at n=33;
   assert the bound, not the measured value).
3. **Sealed limit**: `k_c = 0` reproduces `concentration_bounded()` to 1e-12.
4. **Sink limit**: `Bi = 1e4` gives `C(z_hi)/C(mid) < 1e-3`.
5. **Cross-language anchors** (`Bi=2`, free basis, shipped flow,
   `H=100 µm`, `D_eff=2e-11`, `kH=0.158114`, and
   `Ĉ = 4πD·C/Q`, from a 3000-mode reference):

   | z_s/H | z_t/H | ρ (µm) | Ĉ_robin |
   |---|---|---|---|
   | 0.25 | 0.25 | 2.0 | 5.0484644484e+05 |
   | 0.25 | 0.75 | 10.0 | 7.0837466276e+03 |
   | 0.50 | 0.50 | 20.0 | 5.1114593668e+04 |
   | 0.90 | 0.99 | 5.0 | 1.2272983283e+05 |
   | 0.10 | 0.90 | 40.0 | 4.6022438359e+03 |

6. **CPU/GPU parity** through the shared table, including a `k_c = 0` case so
   the sealed path stays covered on device.

## 9. Approximations retained, to be stated in the docs

- Flow is uniform and evaluated at the source. `radial_velocity` varies with `z`,
  so both the image construction and this eigenmode family are exact only for
  `z`-uniform flow. Unchanged by this work; it is now the largest remaining
  approximation in the kernel and should be quantified separately.
- Correction table keying and construction use the peristalsis-free mean radial
  and distal profiles. This deliberately neglects peristaltic modulation in the
  correction; the sealed image base and runtime drift factor continue to use
  full time-varying flow. The Robin tests report the measured correction change
  at factors 0.5, 1.0, and 1.5 without assigning a pass/fail tolerance. At the
  shipped parameters and a near-lumen target, the maximum measured change was
  `6.00591e-3` relative to the reconstructed field.
- Lateral periodicity is unchanged.
- `δ` is colon-adjusted from small-bowel data, not measured in colon.
