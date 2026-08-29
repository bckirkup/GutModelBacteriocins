# Robin lumen boundary for the QSSA toxin Green's function — design spec

Status: authored by the lead, not yet implemented. Numbers below are measured,
not assumed; the scripts that produced them are named per section.

## 1. What is being replaced and why

The corrected two-wall series from #359 enforces no flux at **both** walls, i.e.
a sealed slab. `AdvectionField::radial_velocity` is nonzero at `z_hi`, so the
lumen wall is not physically sealed, and the choice matters far more than the
bug #359 fixed: measured at shipped parameters (`H=100 µm`, `D_eff=2e-11 m²/s`,
`kH=0.158`), a transparent lumen gives 0.35–0.67× and a perfect sink 0.01–0.52×
the sealed toxin concentration (`lumen_bc_bracket.py`), against a 9–19% bias
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
`b = (k_c + U_z/2)/D`, the two conditions become

```
z_lo:  ψ' + a ψ = 0
z_hi:  ψ' + b ψ = 0
```

so the transformed problem is Robin at **both** walls even though the physical
condition at `z_lo` is pure no-flux. Ignoring this would be a `Pe/2 ≈ 0.043`
error in the mode structure at shipped parameters — small, but it costs nothing
to carry, because the eigenvalue solve happens once at init.

Note `b - a = k_c/D` exactly: the drift cancels from the eigencondition
numerator and survives only through the `ab` product.

## 4. Eigenmodes

```
φ_n(z) = cos(β_n ζ) - (a/β_n) sin(β_n ζ),      ζ = z - z_lo
tan(β_n H) = (b - a) β_n / (β_n² + a b)
N_n = ∫₀^H φ_n(ζ)² dζ
κ_n = sqrt( (λ + |U|²/(4D)) / D + β_n² )

C(r) = exp(U·(x - x_s)/(2D)) · (Q / 2πD) · Σ_n φ_n(ζ_s) φ_n(ζ_t) K₀(κ_n ρ) / N_n
```

`ρ` is the lateral (xy) separation. Since `b - a = k_c/D > 0` and `ab ≥ 0`, the
right-hand side of the eigencondition is positive, so there is exactly one root
per branch `(nπ, nπ + π/2)` including `n = 0`; bracket-and-bisect per branch is
sufficient and cannot miss or double-count a root. As `k_c → 0` the `n = 0` root
goes to zero and the expansion degenerates to the sealed case — handle
`β_0 H < 1e-8` by the limit `φ_0 = 1 - a ζ` rather than dividing by `β_0`.

Validated in `robin_mode_expansion.py`: at `Bi → 0` the expansion reproduces the
sealed image series to 2e-16 relative across nine `(z_s, z_t)` pairs, and
`C(z_hi)/C(H/2)` falls 0.71 → 0.0002 as `Bi` goes 0 → 1e4.

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
of the field only halves the count (199 at ρ=1 µm, `robin_correction_convergence.py`),
so that is not a way out either.

## 6. Implementation: sealed images + tabulated correction

Keep #359's image series for the sealed part and add a per-species table of

```
ΔĈ(z_s, z_t, ρ) := (4πD/Q) · [ C_robin_modes - C_sealed_modes ]
C_robin ≈ C_sealed_images + (Q/4πD) · ΔĈ_interp
```

`ΔC` is regular at `ρ → 0` (the free-space singularity is identical in both
fields and cancels), bounded, and smooth in all three arguments, so it
tabulates on a uniform grid **including** `ρ = 0` — which is exactly what the
full field cannot do. Runtime cost is one trilinear interpolation, identical on
host and device, with no Bessel evaluation in the hot loop.

- Table axes: `z_s, z_t ∈ [z_lo, z_hi]`, `ρ ∈ [0, cutoff]`, uniform, `n = 33`
  per axis (0.27 MB per species in doubles).
- Built once at init from the mode sum with a high mode count (target 2000+;
  `ΔĈ` is converged by 600 at ρ=0.6 µm, verified in `robin_table_accuracy.py`).
- One eigenvalue set per `z_s` row: `a` and `b` depend on `U(z_s)`, and flow is
  already evaluated at the source, so this is consistent with the existing
  approximation rather than a new one.
- Scale by `Q` at runtime; the tabulated quantity is `Q`-independent.

Measured interpolation error against the direct mode sum over 300 random points
(`robin_table_accuracy.py`, `Bi = 2`), as a fraction of the **local total
field**: 4.2e-3 at n=17, 2.8e-3 at n=33, 2.3e-3 at n=65. n=33 is the pick: the
residual is two orders below the 2–3× physics effect this whole change is about,
and the n=65 table costs 8× the memory for no useful gain.

## 7. Parameter: expose the length, not the coefficient

Config key `toxin.lumen_transfer_length` = `δ` [m], with

```
k_c = D_free / δ        (per species, so it scales with each toxin's D_free)
Bi  = k_c·H/D_eff
```

`δ` is the measurable quantity (a lumen-side unstirred/diffusive boundary layer)
and making it the input keeps the two species consistent automatically.

Default `δ = 100 µm` → `Bi ≈ 2`. Basis: human jejunal unstirred-layer thickness
is 35–40 µm (Levitt et al., *J Clin Invest* 1990, range 23–65; Strocchi &
Levitt, *Am J Physiol* 1992, ≈35 µm), which replaced older ~600 µm osmotic
estimates. Those are small-bowel measurements with glucose probes; the colon is
less stirred and carries an outer mucus layer, so the default is deliberately
2–3× the jejunal value and is documented as **colon-adjusted, not measured**.

`Bi = 0.3 / 2 / 5` (δ ≈ 600 / 100 / 40 µm) puts near-lumen toxin at ≈0.75× /
0.45× / 0.2× the sealed value, so `δ` must appear in the sensitivity sweep for
any lethal-core/halo or comet-tail claim regardless of the default.

Limits to keep reachable and tested: `δ → ∞` (`k_c = 0`) must return the #359
sealed result bit-for-bit via the same code path, and large `Bi` must drive
`C(z_hi) → 0`.

## 8. Required tests (independent oracles, not self-oracles)

1. **Flux residual at both walls**, numerically differentiated: `-D ∂C/∂z + U_z C = 0`
   at `z_lo` and `-D ∂C/∂z - k_c C = 0` at `z_hi`, over several source
   positions and lateral offsets, normalised by `D·C/H`.
2. **Table vs direct mode sum**: max relative error over random points below
   5e-3 of the local total field (measured 2.8e-3 at n=33; assert the bound, not
   the measured value).
3. **Sealed limit**: `k_c = 0` reproduces `concentration_bounded()` to 1e-12.
4. **Sink limit**: `Bi = 1e4` gives `C(z_hi)/C(mid) < 1e-3`.
5. **Cross-language anchors** (`Bi=2`, `H=100 µm`, `D_eff=2e-11`, `kH=0.158114`,
   `ΔĈ` and `Ĉ = 4πD·C/Q`, from a 3000-mode reference):

   | z_s/H | z_t/H | ρ (µm) | Ĉ_robin | ΔĈ |
   |---|---|---|---|---|
   | 0.25 | 0.25 | 2.0 | 5.1061945263e+05 | -3.3149597263e+04 |
   | 0.25 | 0.75 | 10.0 | 1.9964906375e+04 | -3.7235028089e+04 |
   | 0.50 | 0.50 | 20.0 | 5.0272543392e+04 | -3.6163301346e+04 |
   | 0.90 | 0.99 | 5.0 | 1.2537873736e+05 | -7.7533465488e+04 |
   | 0.10 | 0.90 | 40.0 | 1.0721327524e+04 | -3.7660834321e+04 |

6. **CPU/GPU parity** through the shared table, including a `k_c = 0` case so
   the sealed path stays covered on device.

## 9. Approximations retained, to be stated in the docs

- Flow is uniform and evaluated at the source. `radial_velocity` varies with `z`,
  so both the image construction and this eigenmode family are exact only for
  `z`-uniform flow. Unchanged by this work; it is now the largest remaining
  approximation in the kernel and should be quantified separately.
- Lateral periodicity is unchanged.
- `δ` is colon-adjusted from small-bowel data, not measured in colon.
