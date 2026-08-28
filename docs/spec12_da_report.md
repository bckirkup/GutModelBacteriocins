# Spec 12 reaction–diffusion scale-gating report

## Definitions and numerical method

I solved the dimensionless spherical Michaelis–Menten boundary-value problem with SciPy `solve_bvp` and exact matching to the infinite exterior. The exterior matching condition is `c'(1)=1-c(1)`. The effectiveness factor is normalized to the uptake that would occur if every cell experienced the bulk concentration. Michaelis–Menten kinetics has no finite mathematical anoxic front: near zero concentration it becomes first order, so concentration approaches zero exponentially. I therefore define an operational dead zone by `C/Cinf <= 1e-3`. The planar penetration depth uses the same threshold.

## Baseline spherical results

Oxygen baseline: D=2.0e-9 m2/s, Cinf=20 µM, Km=0.5 µM, q=8.33e-19 mol/cell/s. Carbon baseline: D=6.0e-10 m2/s, Cinf=100 µM, Km=sqrt(10)=3.16 µM (geometric midpoint of 1–10 µM), and recomputed q=8.33e-19 mol/cell/s.

### oxygen

| density (cells/mL) | Da=1 (µm) | R10 (µm) | R50 (µm) | R90 (µm) | R(C/Cinf=1e-3) (µm) | planar 1e-3 depth (µm) |
|---:|---:|---:|---:|---:|---:|---:|
| 1e+09 | 219.13 | 99.33 | 223.81 | 310.63 | 392.72 | 430.87 |
| 1e+10 | 69.30 | 31.41 | 70.78 | 98.23 | 124.19 | 136.25 |
| 1e+11 | 21.91 | 9.93 | 22.38 | 31.06 | 39.27 | 43.09 |
| 1e+12 | 6.93 | 3.14 | 7.08 | 9.82 | 12.42 | 13.63 |

### carbon_glucose

| density (cells/mL) | Da=1 (µm) | R10 (µm) | R50 (µm) | R90 (µm) | R(C/Cinf=1e-3) (µm) | planar 1e-3 depth (µm) |
|---:|---:|---:|---:|---:|---:|---:|
| 1e+09 | 268.33 | 122.05 | 275.56 | 385.66 | 506.86 | 554.83 |
| 1e+10 | 84.85 | 38.60 | 87.14 | 121.96 | 160.28 | 175.45 |
| 1e+11 | 26.83 | 12.21 | 27.56 | 38.57 | 50.69 | 55.48 |
| 1e+12 | 8.49 | 3.86 | 8.71 | 12.20 | 16.03 | 17.55 |

## Sensitivity ranges

Ranges below are full factorial extrema, not confidence intervals. They deliberately combine all requested extremes.

| medium | density | R50 min–max (µm) | R(dead) min–max (µm) |
|---|---:|---:|---:|
| oxygen | 1e+09 | 24.92–1266.60 | 42.06–3091.66 |
| oxygen | 1e+10 | 7.88–400.53 | 13.30–977.67 |
| oxygen | 1e+11 | 2.49–126.66 | 4.21–309.17 |
| oxygen | 1e+12 | 0.79–40.05 | 1.33–97.77 |
| carbon_glucose | 1e+09 | 29.06–2706.34 | 73.05–4194.76 |
| carbon_glucose | 1e+10 | 9.19–855.82 | 23.10–1326.50 |
| carbon_glucose | 1e+11 | 2.91–270.63 | 7.31–419.48 |
| carbon_glucose | 1e+12 | 0.92–85.58 | 2.31–132.65 |

The cell-specific uptake rate dominates the oxygen radius uncertainty among the requested inputs. One factor at a time around the baseline, the R50 spans are 10.0-fold for q, 3.09-fold for Cinf, 1.45-fold for D, and 1.19-fold for Km (these ratios do not depend on density). For carbon, the requested Cinf and q sweeps give 8.06-fold and 10.0-fold R50 spans, respectively, while Km gives 1.07-fold; thus q and the poorly constrained carbon concentration are comparably dominant. The larger tabulated full-factorial envelopes combine all extremes and are not confidence intervals.

## Carbon-rate correction

`10 mmol glucose gDW^-1 h^-1 × 3e-13 gDW cell^-1 = 3e-15 mol cell^-1 h^-1 = 8.33e-19 mol cell^-1 s^-1`, not `2.8e-19`. I used `8.33e-19` as the carbon baseline and factorial center. The supplied `2.8e-19` value is retained as separately tagged rows in the CSV. Carbon concentration is explicitly treated as poorly constrained.

## Closed-form comparisons

### zero-order test (k=1e-4)

| depletion | MM A | zero-order radius error | first-order radius error |
|---:|---:|---:|---:|
| 0.1 | 0.2 | -0.0% | -99.0% |
| 0.5 | 1 | -0.0% | -98.7% |
| 0.9 | 1.799 | +0.0% | -97.8% |

### baseline oxygen (k=0.025)

| depletion | MM A | zero-order radius error | first-order radius error |
|---:|---:|---:|---:|
| 0.1 | 0.2055 | -1.3% | -83.7% |
| 0.5 | 1.043 | -2.1% | -79.6% |
| 0.9 | 2.009 | -5.4% | -66.6% |

### first-order test (k=10)

| depletion | MM A | zero-order radius error | first-order radius error |
|---:|---:|---:|---:|
| 0.1 | 2.382 | -71.0% | -4.3% |
| 0.5 | 18.34 | -76.6% | -2.8% |
| 0.9 | 91.29 | -86.0% | -0.9% |

At `Da=1`, the numerical baseline oxygen result is 47.99% center depletion with eta=0.9848; carbon is 47.50% with eta=0.9813. The zero-order spherical solution before core exhaustion is `C(0)/Cinf = 1 - Da/2`; therefore `Da=1` is exactly the 50% center-depletion radius in the zero-order model, not the onset of depletion. In the low-concentration first-order limit, `C(0)/Cinf = 1/cosh(sqrt(Da/k))`. The numerical Michaelis–Menten solution converges to each expression in its respective limit.

**Spec rule of thumb:** `R50 ≈ sqrt(D*Cinf/(rho*q_cell))` when `Km << Cinf`; interpret this as the radius for ~50% center depletion, not an anoxic-core threshold. Use the Michaelis–Menten multiplier `sqrt(A_target(Km/Cinf))` for other depletion targets; a practical dead-zone threshold is `C(0)/Cinf <= 1e-3`.

## Numerical convergence

The worst relative difference between the two finest runs across the three verification cases was 1.96e-10. Full results are in `spec12_da_convergence.csv`.

## Assumptions and limitations

- All kinetic and transport values are user-supplied assumptions; I did not fetch literature.
- Uniform cell density, constant diffusivity, steady state, spherical symmetry, and an infinite external medium are assumed.
- The slab result is semi-infinite and one-sided with a fixed source concentration. A finite mucus slab or a second boundary changes the profile.
- “Dead zone” is threshold-dependent because Michaelis–Menten uptake does not create a finite exact zero-concentration region.