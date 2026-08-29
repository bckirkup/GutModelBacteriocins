# Design Review: Tabulated Robin Correction for GutIBM Toxin Green Functions

**To:** Benjamin Kirkup and Devin  
**Review basis:** GutModelBacteriocins PR [#359](https://github.com/bckirkup/GutModelBacteriocins/pull/359), merged as [`49dae6827458cca61ee360419c637a6ff17d50b2`](https://github.com/bckirkup/GutModelBacteriocins/commit/49dae6827458cca61ee360419c637a6ff17d50b2)  
**Status:** Conditional approval for a staged prototype

## Executive recommendation

Keep the decomposition

\[
G_{\mathrm{Robin}}=G_{\mathrm{sealed}}+\Delta G,
\qquad
\Delta G=G_{\mathrm{Robin}}-G_{\mathrm{sealed}},
\]

and precompute the correction. This is a sound way to retain the corrected two-wall image series from PR #359 while removing the Robin mode sum from the hot loop. Expose a transfer length \(\delta\), rather than making an opaque transfer coefficient \(k_c\) the ordinary user parameter. However, do not implement the proposal as one three-dimensional table per receptor-target field indexed only by \((z_s,z_t,\log\rho)\). GutIBM transport parameters belong to each bacteriocin locus and source, and the current kernel is advective and direction-dependent. Approval should therefore be conditional on an explicit declaration of what concentration the kernel transports and which diffusivity belongs in the boundary law.

> **Decision**
>
> - **YES:** retain the sealed kernel and add a tabulated correction.
> - **YES:** expose transfer length \(\delta\) in ordinary configuration.
> - **NO:** do not use a simple per-field-species \((z_s,z_t,\log\rho)\) table.
> - **CONDITIONAL:** proceed only after fixing the concentration-state definition, boundary diffusivity convention, advection treatment, and table key.

## 1. Why the subtraction is attractive

For a source and target strictly inside the slab, the Robin and sealed Green functions have the same local free-space term, \(Q/(4\pi D_{\mathrm{eff}}r)\). Their difference cancels the \(1/r\) singularity. At fixed positive distance from each boundary and with a screened steady operator, \(\Delta G\) is bounded and smooth through \(\rho=0\), where \(\rho\) is lateral source-target separation. That makes the correction substantially easier to interpolate than either full kernel.

This also addresses the immediate cost problem. The corrected implementation already has a shared host/device two-family Neumann sum, a relative tolerance of \(10^{-10}\), and an adaptive cap of 512 shells ([`src/diffusion/neumann_image_series.h:10-20,58-107`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/diffusion/neumann_image_series.h#L10-L107)). In the difficult short-range case cited in the proposal, a direct Robin evaluation can require about 199 mode terms at \(\rho=1\,\mu\mathrm{m}\). A table lookup replaces that variable work with a fixed number of loads and fused interpolation operations, with no Bessel evaluation in the runtime loop.

Use one table representation and one interpolation specification on host and device. The requirement should be the **same data, interpolation contract, and numerical tolerance**, not bitwise identity. The prior CUDA audit found genuine device kernels but also floating-point atomic accumulation and architecture-dependent ordering, so universal bitwise parity is not an established property of this codebase (task:180ea40c-de07-4423-ada9-8675f07c6e84). Preserve the direct mode sum as an independent oracle. It should not call the table or share interpolation code.

PR #359 is an appropriate sealed baseline. It replaced the earlier incomplete image construction identified by the scientific audit (task:cdd046b0-61fe-459b-9725-e2b2364d6654) with both translated and reflected source families and records convergence-cap hits. Do not duplicate or fork that series inside the Robin feature.

## 2. Scope of the regularity claim

The useful regularity statement needs two qualifications.

### Near the Robin wall

If both source and target approach the Robin boundary while \(\rho\to0\), the high-mode difference between Robin and Neumann reflection coefficients decays too slowly to guarantee a bounded axis limit. A logarithmic near-axis singularity can remain in \(\Delta G\). This is not merely a mathematical corner case because agents may be clamped to a wall.

The implementation therefore needs one of the following, selected and tested explicitly:

1. enforce a minimum source and target interior distance tied to cell radius or grid geometry;
2. use a direct Robin branch inside a near-wall/near-axis guard region; or
3. tabulate a separately regularized asymptotic remainder after subtracting the boundary singular term.

Option 2 is the safest prototype.

### Vanishing screening

For a steady operator with negligible decay and negligible effective advective screening, the sealed slab has a non-removing zero mode. Its steady Green function diverges. A Robin lumen can remain finite because it removes mass. In that limit, writing a finite Robin result as sealed plus a negative, divergent correction creates catastrophic cancellation, or the two terms are individually undefined.

Add a direct Robin path or reject this regime based on a documented threshold in \(\kappa H\) and the relevant advective screening measure. The corrected image code already treats scaled screening below 0.05 as a forced-cap case ([`src/diffusion/neumann_image_series.h:27-52`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/diffusion/neumann_image_series.h#L27-L52)); the Robin decomposition needs a stronger accuracy-based gate because cancellation, not only series truncation, is at issue.

**Therefore, “\(\Delta G\) is regular” is an interior, screened-operator result. It is not a global theorem over the closed slab.**

## 3. Table keys must follow transport classes, not fields

GutIBM stores free diffusivity, retardation, and protease half-life in each `BICluster` ([`src/core/types.h:95-113`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/core/types.h#L95-L113)). The quasi-steady-state approximation (QSSA) solver copies those locus values into each source's `GreensFunctionParams` and adds source-position-dependent dilution to degradation ([`src/diffusion/qssa_solver.cpp:322-351`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/diffusion/qssa_solver.cpp#L322-L351)). It then groups sources by receptor target when depositing receptor-targeted fields ([`src/diffusion/qssa_solver.cpp:641-653`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/diffusion/qssa_solver.cpp#L641-L653)).

This proves that field species and transport class are different concepts. ColE1 and ColE2 both target BtuB, but the plasmid library gives them different \(D_{\mathrm{free}}\), retardation, and toxin properties ([`src/genome/plasmid.cpp:37-70`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/genome/plasmid.cpp#L37-L70)). One BtuB correction table would apply the wrong kernel to at least one source.

Normalize coordinates by mucus height \(H\), and cache tables by a nondimensional transport tuple such as

\[
\left(\mathrm{Bi},\ \kappa H,\ \mathrm{Pe}_{\parallel},\
\mathrm{Pe}_{n},\ \text{boundary-mode pair}\right),
\]

with geometry \(H\) and the normalization convention recorded. For constant local coefficients,

\[
\kappa^2=\frac{\lambda}{D_{\mathrm{eff}}}
          +\frac{|\mathbf U|^2}{4D_{\mathrm{eff}}^2},
\quad
\mathrm{Pe}_{\parallel}=\frac{U_{\parallel}H}{D_{\mathrm{eff}}},
\quad
\mathrm{Pe}_{n}=\frac{U_nH}{D_{\mathrm{eff}}}.
\]

Here \(\kappa\) is the screening parameter of the gauge-transformed Helmholtz operator, and the signed Peclet numbers encode lateral and normal drift. “Boundary-mode pair” distinguishes, for example, sealed epithelium plus Robin lumen from other future combinations.

If the plasmid library remains a small exact set, constructing one table for each distinct source parameter tuple is simpler and less error-prone than a high-dimensional universal table. If configuration overrides or mutations can vary transport continuously, use a quantized cache with a declared interpolation/error policy and a maximum cache size. Do not create an unbounded new table for every floating-point tuple.

## 4. Advection cannot be collapsed into scalar \(\rho\)

The current free kernel contains

\[
\exp\!\left[\frac{\mathbf U\cdot(\mathbf r_t-\mathbf r_s)}
{2D_{\mathrm{eff}}}\right]
\]

along with isotropic screening. In code, velocity is sampled at the source, and the signed dot product sets downstream versus upstream concentration ([`src/diffusion/greens_function.cpp:275-300,303-317`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/diffusion/greens_function.cpp#L275-L317)). A correction indexed only by scalar \(\rho\) discards that signed direction.

Prefer a gauge transform. For locally constant \(\mathbf U\), write

\[
C(\mathbf r)=
\exp\!\left[\frac{\mathbf U\cdot(\mathbf r-\mathbf r_s)}{2D_{\mathrm{eff}}}\right]
\psi(\mathbf r),
\]

then table the smoother Helmholtz correction for \(\psi\) and restore the exponential at lookup. The transformed Robin coefficient must include normal velocity. With outward-normal convention \(-D_{\mathrm{eff}}\partial_n C=k_c C\), the transformed condition is

\[
-D_{\mathrm{eff}}\partial_n\psi
=\left(k_c+\frac{U_n}{2}\right)\psi.
\]

The sign must be re-derived if the code uses the opposite normal or flux convention.

Advection in GutIBM varies with depth and time. The profile depends on \(z\), and the peristaltic factor depends on current time and, when wavelength is positive, axial position ([`src/fields/advection.cpp:12-18,31-50`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/fields/advection.cpp#L12-L50)). A table built once from \(\mathbf U(z_s,t=0)\) becomes stale.

Rank the implementation options as follows:

1. **A. Universal nondimensional table interpolated in Peclet number.** Recommended if error maps show a tractable sparse or tensor representation.
2. **B. Small velocity-bin cache.** Acceptable when actual peristaltic extrema can be covered by a bounded set of bins with verified interpolation error.
3. **C. Diffusion-only correction.** Acceptable only with an explicit, validated error envelope over configured flow and source depths.

The existing source-local constant-velocity approximation is itself a documented scientific risk (task:cdd046b0-61fe-459b-9725-e2b2364d6654). This feature should not silently add a second approximation by dropping directionality.

## 5. Diffusivity, transfer length, and the Biot number

The current kernel computes

\[
D_{\mathrm{eff}}=\frac{D_{\mathrm{free}}}{R}
\]

for each source ([`src/diffusion/greens_function.cpp:303-317`](https://github.com/bckirkup/GutModelBacteriocins/blob/49dae6827458cca61ee360419c637a6ff17d50b2/src/diffusion/greens_function.cpp#L303-L317)). If the boundary law is

\[
-D_{\mathrm{eff}}\partial_n C=k_c C,
\qquad k_c=\frac{D_{\mathrm{free}}}{\delta},
\]

then

\[
\mathrm{Bi}=\frac{k_cH}{D_{\mathrm{eff}}}
=\frac{RH}{\delta},
\]

not \(H/\delta\). For \(H/\delta=2\), the supplied defaults produce:

| Toxin | Retardation \(R\) | Realized Bi |
|---|---:|---:|
| ColB | 1.27 | 2.5 |
| Microcin V | 1.23 | 2.5 |
| ColE2 | 2.04 | 4.1 |
| ColIa | 5.17 | 10.3 |
| ColE1 | 50.22 | 100.4 |
| ColM | 55.15 | 110.3 |

These are values supplied in the proposal pending citation and configuration verification. They show that “Bi \(\sim2\)” applies only to weakly retarded toxins. Under this convention, the lethal-core toxins approach an absorbing lumen.

Two model choices are coherent, but they represent different state variables:

- **Choice F, free concentration.** \(C\) is free toxin. An external-film relation \(k_c=D_{\mathrm{free}}/\delta\) can be appropriate, but retardation or partitioning must be represented explicitly, and the outlet flux must act on the free fraction.
- **Choice T, retarded effective concentration.** \(C\) is the total or effective transported concentration. Use \(k_{c,\mathrm{eff}}=D_{\mathrm{eff}}/\delta\), or derive a partition-adjusted coefficient from the chosen sorption model. Then \(\mathrm{Bi}=H/\delta\) across species.

Write a one-page state-variable declaration before coding. It should define the units and physical meaning of \(C\), source \(Q\), retardation, decay, free fraction, boundary flux, and receptor exposure. Do not choose the convention because it makes table reuse convenient.

## 6. Evidence for \(\delta\)

The proposal supplies human jejunal unstirred-layer glucose estimates of 35 to 40 µm, with a stated range of 23 to 65 µm. These are **values supplied in the proposal pending citation verification**. They concern small-bowel, small-solute transport, not colonic bacteriocin proteins. They do not establish a colon toxin transfer length.

A default \(\delta=100\,\mu\mathrm{m}\) is acceptable as an explicit colon-adjusted prior for exploration, not as an empirical colon default. Sweep broadly, for example 30, 100, and 300 µm, or use a dimensionless Bi grid. Every run should output realized Bi for each transport tuple.

Likewise, attenuation ratios 0.75, 0.45, and 0.2 at Bi 0.3, 2, and 5 are geometry-specific. Whenever quoted, include source depth, target depth, \(\rho\), \(\kappa H\), flow, and the concentration normalization. Exact citations should be verified before final parameter documentation unless actual BibTeX is available. No bibliography should be inferred from surnames or secondary summaries.

## 7. Table and interpolation design

Use the following initial design, then adapt from measured error maps:

- Include an explicit \(\rho=0\) plane.
- Near the axis, interpolate against \(s=\rho^2\), which respects even radial smoothness for interior points.
- Beyond a matched radius, switch to \(\log\rho\). Match value and preferably first derivative at the transition.
- Exploit \(z_s,z_t\) reciprocity only where the transformed operator and boundary conditions support it. With advection, store or enforce the correct gauge-weighted reciprocity, not naive symmetry.
- Validate the reconstructed total \(C_{\mathrm{Robin}}\), not only \(\Delta C\). Subtraction can have a small correction error yet a large total relative error near cancellation.
- Preserve nonnegative total concentration. Do not simply clip as a substitute for accurate interpolation; positivity failures should be diagnostics and acceptance failures.
- Check the lumen residual \(D\partial_n C+k_cC\) using the exact code sign convention, and the epithelial Neumann residual.
- Record table version, grid, nondimensional key, mode-sum tolerance, interpolation method, state-variable basis, and a cryptographic hash in run provenance.

Start with \(N_s=N_t=48\) to 64 and \(N_\rho=128\) to 256. Use adaptive refinement in the wall, axis, and cancellation regions before compression.

| Grid | Double-precision values | Memory per parameter tuple |
|---|---:|---:|
| \(48\times48\times128\) | 294,912 | 2.25 MiB |
| \(48\times48\times256\) | 589,824 | 4.50 MiB |
| \(64\times64\times192\) | 786,432 | 6.00 MiB |
| \(64\times64\times256\) | 1,048,576 | 8.00 MiB |

The 8 MiB figure is per **parameter tuple**, not necessarily per field species. A naive Cartesian extension across Bi, \(\kappa H\), and two Peclet dimensions would grow rapidly. Establish required ranges from real configurations first, then consider low-rank compression, sparse grids, or bounded caches only after error mapping.

## 8. Concrete implementation decision tree

1. **Declare the state variable.** Specify free versus total/effective toxin and the free fraction seen by receptors and the lumen.
2. **Select the boundary diffusivity convention.** Derive \(k_c\), Bi, and flux units consistently from Step 1.
3. **Derive the transformed mode sum.** Include decay, signed parallel flow, normal velocity in the transformed Robin coefficient, and the exact normal convention.
4. **Define keys and ranges from actual inputs.** Inventory the plasmid library, configuration overrides, mutation paths, \(H\), decay, and peristaltic extrema.
5. **Build the direct oracle and table.** Keep the direct Robin mode sum independently callable. Build tables offline or once at initialization, with deterministic metadata and a hash.
6. **Implement host interpolation first.** After it passes, copy the same table and contract to CUDA.
7. **Gate unsupported regimes.** Guard near-wall/near-axis points, low screening, out-of-range keys, excessive cancellation, and stale velocity bins with direct or reject policies.
8. **Benchmark after accuracy.** Report lookup latency, cache behavior, memory, and end-to-end kernel time. Preserve acceptance tolerances.

## 9. Acceptance test matrix

| Test axis | Required cases and acceptance evidence |
|---|---|
| Independent oracle | Random off-grid comparisons with the direct mode sum across \(z_s,z_t,\rho\); table and oracle must not share interpolation code. |
| Axis | Exact \(\rho=0\), both sides of the matched radius, and logarithmically spaced near-axis points. |
| Walls | Source and target near the Robin lumen, near the sealed epithelium, and interior; include the direct-branch threshold. |
| Boundary strength | Bi = 0, 0.3, 2, 5, 20, 100, plus the Bi \(\to\infty\) Dirichlet limit. Bi = 0 must recover sealed behavior. |
| Screening | Low, normal, and high \(\kappa H\); verify direct/reject behavior in the ill-conditioned low-screening regime. |
| Flow | Positive and negative parallel Peclet number, nonzero normal Peclet number, zero flow, and peristaltic velocity extrema. |
| Transport identity | ColE1 and ColE2 deposited to the same BtuB-targeted field but evaluated with distinct source transport tuples. |
| Boundary residuals | Robin residual at the lumen and Neumann residual at the epithelium, evaluated away from the point source. |
| Physical totals | Total positivity, absolute and relative interpolation error, cancellation condition number, and error in total Robin kernel, not only correction. |
| CPU/CUDA | Same table hash and parameter metadata; tolerance parity over identical query points. Do not require universal bitwise equality. |
| Performance | Confirm fixed-cost \(O(1)\) lookup, report host and device latency and memory, and show no relaxation of accuracy criteria. |

Error thresholds should combine relative and absolute tolerances because the total may approach zero. Report percentile and worst-case errors by wall distance, axis distance, screening, and flow bin rather than only one global maximum.

## 10. Draft configuration surface

The following names are a draft, not an approved schema:

```ini
toxin_boundary.mode = "robin"
toxin_boundary.transfer_length = 1.0e-4
toxin_boundary.concentration_basis = "free" | "effective"
toxin_boundary.table.enabled = true
toxin_boundary.table.rho_points = 192
toxin_boundary.table.z_points = 64
toxin_boundary.table.rel_tolerance = 1e-5
toxin_boundary.table.low_screening_policy = "direct" | "reject"
```

Do not expose \(k_c\) independently in the ordinary interface. If an expert override is retained, record both \(k_c\) and \(\delta\), require an explicit basis/diffusivity convention, and reject simultaneous inconsistent settings rather than selecting one silently. Also record the realized Bi values after all per-locus parameters and overrides are resolved.

## Final recommendation to Devin

Implement this in two pull requests.

**PR A: mathematics, oracle, and CPU table.** Turn advection and time-varying flow off. Add the direct Robin mode sum, low-screening and near-wall guards, a CPU correction table, provenance, and the full diffusion-only boundary/interpolation test matrix. Demonstrate accuracy of the total kernel and both boundary residuals before optimizing.

**PR B: advection, CUDA, and peristaltic parameterization.** Add the gauge-transformed advective correction, Peclet parameterization or a validated velocity-bin cache, device lookup using the exact table hash, and CPU/CUDA tolerance tests at peristaltic extrema.

Do not combine the table mathematics and GPU implementation in one pull request. Separating them isolates analytical and interpolation correctness from acceleration, makes review failures attributable, and prevents a fast implementation from becoming the de facto oracle.
