# External audit of 2026-08-29 — claims, verdicts, and what we decided

An external automated review (KOSMOS, Edison Scientific platform) audited the
repository at `7bfa7c6`. This file is the record: what was claimed, which claims
survived being checked against source, which did not, and what we decided as a
result. It exists so the reasoning outlives the conversation it happened in, and
so a future reader can tell a verified defect from an accepted interpretation.

Verdicts were established against `main` @ `89c40bf`; `7bfa7c6` is its parent and
no executable source differs between them.

## Standing of the report as a whole

Of the 16 ranked defects checkable against source, 14 were factually correct as
stated, one was correct in conclusion but wrong in mechanism, and one was correct
but understated. No file or line citation was fabricated. The report's weakness
was not accuracy but novelty and framing: several findings were already recorded
here as known open decisions, and one recommendation collided with a deliberate
guard already in the parser.

Treat the report as a competent code review, not as an authority on this model's
scientific scope. Nothing in it required a third opinion; what remains open
requires experiments, not arbitration.

## Verified defects

| # | Claim | Verdict |
|---|---|---|
| 1 | Two-wall Neumann image series incomplete and duplicated | Confirmed, measured, and gated — see below |
| 2 | `uptake_limit` ships as `none`, so growth is unfunded by realized removal | Confirmed factually; the remedy is a decision, not a patch — see below |
| 3 | Configuration fails open (missing file, malformed JSON, unknown keys) | **Confirmed and fixed.** JSON-vs-legacy selection now uses file shape rather than an error-message whitelist: once a root JSON object is found, every parse failure is fatal; files without one retain legacy flat-key parsing. Named-file open failures, malformed arrays, unknown keys, and invalid known values fail closed by default. `GUTIBM_STRICT_CONFIG` is inverted: unset/`1` is strict, while `0` is the explicit lenient escape hatch. |
| 4 | GPU per-species diffusion completion suppresses CPU fallback | **Confirmed and fixed.** A whole-field pre-flight now declines GPU diffusion before any launch when an eligible species exceeds the line cap, so the pipeline's existing host fallback covers every species. The mixed-eligibility window is exactly `nz == 1025`: Dirichlet solves `nz - 1`, Robin/Flux solve `nz`; above the cap every species is refused uniformly. Their `nz = 1025` reproducer is covered host-side for the predicates and on the T4 for the device path; a species the host does not diffuse never forces a fallback |
| 5 | GPU mechanics omits CDI corpses | **Confirmed and fixed.** Host and device now share the same death-time/persistence predicate; fresh CDI corpses remain mechanical obstacles while expired corpses are excluded. Host-side coverage is local, while the device assertion is evidenced only by the T4 CI job |
| 6 | Checkpoint restart does not restore RNG state | Confirmed, low novelty — `docs/BRANCHING_FROM_CHECKPOINTS.md` already states a fork is a population-state continuation, not a bit-identical restart |
| 7 | Post-division SOS hazard applies to both parent and daughter | Confirmed: 1.98% per division event against a documented 1% "per division". An estimand ambiguity, cheap to settle |
| 8 | Comet-tail analysis fabricates geometry | **Confirmed and fixed.** Validation now uses physical HDF5 cell centers in `(z, y, x)` flattening order, the model's positive-x distal-flow convention, and an explicit producer reference with periodic minimum-image wrapping |
| 9 | Resident-retention estimator measures lineage-label persistence | Confirmed |
| 10 | Flagship diversity scenario has no immigration block and no grid output | Confirmed, both halves, against a README that describes periodic lumen immigration |
| 11 | Taylor–Aris dispersion toggle is unused | **Confirmed and resolved by removal.** The dead, default-on toggle and uncalled implementation were removed rather than wired into the isotropic transport model. At `z=h`, the measured `D_taylor/D` enhancements are `6.379e-9` for oxygen (`D=2e-9`), `6.379e-5` for carbon (`D=2e-11`), `4.023e-2` for shipped ColE1-like colicin (`D=7.964e-13`), and `2.552` for a hypothetical `D=1e-13`; the enhancement falls as `(z/h)^3`. The long-time limit is valid for the shipped colicin (`t/t_diff=3.44`) but not for the large-effect hypothetical case (`0.43`). The `210` prefactor is also the Poiseuille result while the shipped profile exponent is `1.5`. Honest wiring would require an anisotropic streamwise transport kernel and re-derivation of both the image series and Robin correction table. Reinstating it requires that work. |
| 12 | VBF drag and carrying capacity unused | Half wrong, right conclusion — see below |
| 13 | Requested HDF5 output can fail open | **Confirmed and fixed.** Requested HDF5 output now throws before compute when path validation or file creation fails; under MPI, rank 0 makes the decision and broadcasts it so every rank fails together rather than continuing without the requested record. The failed output file is no longer silently removed, and the error names the file plus the underlying reason where available. |
| 14 | Provenance and compiler contract incomplete | Confirmed: git SHA falls back to `unknown-git-unavailable`, and `<format>` is used widely while the README never states the GCC 13+ requirement |
| 15 | Python manifest writes are non-atomic | Confirmed: destination opened in `"w"`, no temp file, no fsync, no retained generation |
| 16 | Ethanolamine absolute units off by 1000 | Confirmed, already recorded in `docs/UNITS_AUDIT.md`. `eut_km` is off by the same factor, so the Monod penalty is numerically unchanged — labels and coupling only |

### Claim 8 measurements and resolution

The shipped `examples/eari_vadi_validation` run at `step_000015` was measured
with both the former fabricated line and physical HDF5 cell centers:

| Geometry and reference | `comet_tail_ratio` | `comet_tail_asymmetry` |
|---|---:|---:|
| Fabricated geometry | 0.9529 | 0.9927 |
| Real geometry, domain-midpoint reference | 0.7071 | 0.6040 |
| Real geometry, producer reference with minimum-image wrapping | 0.6464 | 0.5923 |

The synthetic oracle makes the geometry defect especially clear: the same
downstream-elongated plume scores `1.0000` / `1.0005` with fabricated
coordinates, but `5.9998` / `35.62` with real coordinates.  Validation now
uses voxel centers derived from HDF5 domain metadata, the model's positive-x
distal-flow convention from `AdvectionField::velocity`, and a documented
producer-reference fallback chain with the x-domain period for wrapping.

With real geometry that shipped example shows **no** comet tail (ratio < 1),
so the `min: 1.5` target has never been met by a real-geometry measurement and
remains unvalidated — we are recording that, not changing the target.

## Where we corrected the report

**#12 VBF drag — right conclusion, wrong mechanism.** `drag_force()` *is* called
from the physics module, so "stored but unused" is inaccurate. The real defect is
different and slightly worse: `a.v` is written only by that drag term, drag is
`-c·v`, and mechanics is overdamped and never writes `a.v`, so a velocity that
starts at zero stays identically zero and the `a.x += a.v·dt` update is a
permanent no-op. `local_capacity()` genuinely has no callers, so the
carrying-capacity half stands as written.

**Daughter placement near bounded walls (Section 4).** `Domain::apply_pbc`
clamps non-periodic axes, so near-wall daughters are never lost from the domain —
they accumulate on the wall. The artifact is a surface density bias, not
immediate loss, which makes the risk milder and the test cheaper than claimed.

## What we found that the report did not

**The CUDA path carried the identical image-series defect.** The same four
duplicated reflected families and the same `N_IMAGES = 3` appeared in
`src/gpu/gpu_common.cuh`. Because both backends were wrong in the same way, no
CPU/GPU parity test could ever have detected it. The correction is now a single
shared host/device series header for exactly that reason.

**Quantified consequence of the image-series defect.** Independent recomputation
reproduced the reported wall-gradient residuals (median 0.044, max 0.216 for a
quantity that should be zero) and added what the report did not measure:
concentration was 9–19% low at `z_s = 0.3H`, and the corrected series needs on
the order of 30 image pairs to converge at default parameters, against the three
that were used.

**Truncation cannot be a fixed constant.** The screening rate sets the shell
count. With decay and flow both negligible the series diverges logarithmically —
a sealed slab with no removal has no bounded steady state, so no shell count is
correct and the configuration is outside the validated envelope. The shipped
implementation therefore derives its shell count from the screening rate, flags
saturation into run provenance, and the parser warns (and refuses under
`GUTIBM_STRICT_CONFIG`) when `sqrt(decay_rate/D_eff)·H < 0.05`. That threshold
deliberately ignores advective screening: `distal_velocity` and
`radial_velocity` both vanish at `z_lo`, so a source at the epithelium is
unscreened by flow no matter how fast the lumen moves.

**The lumen wall should not be reflecting at all.** `radial_velocity` is nonzero
at `z_hi`, so there is advective throughflow at the lumen boundary and
reflecting toxin there contradicts the transport the rest of the model uses. The
epithelium at `z_lo` is a genuine no-flux wall. The corrected two-wall Neumann
series is the right baseline and a prerequisite either way, but the lumen
boundary is being changed to outflow as a separate, subsequent change. Until
that lands, near-lumen toxin footprints are over-retained.

**The `delivery` recommendation collides with a deliberate guard.** The parser
refuses `uptake_limit="delivery"` with `gpu_enabled=true` because CUDA delivery
parity does not exist (`docs/CUDA_DELIVERY_PARITY.md` names the three missing
device pieces). "Make funded uptake mandatory for scientific configurations" is
therefore currently equivalent to "no GPU campaigns" — and the report never
mentions `sherwood`, the unblocked middle path.

## Decisions taken

**Uptake mode.** `sherwood` is the operating mode: `src/fixes/uptake_limit.h` is
a shared `__host__ __device__` closed form and the agent kernel takes
`uptake_limit_mode` straight through, so `none`/`sherwood`/`voxel` all run on
device and only `delivery` falls back to host. Two qualifications must travel
with that decision rather than being discovered later:

- `sherwood` caps demand at `4πD·r·C` and scales `mu_realized` by the same
  fraction, which removes the unbounded-growth artifact. It does **not** deliver
  mass closure: removal remains a prescribed voxel reaction that can be clipped,
  and `closure.enforce_reaction_clip` ships `false`. Growth can still exceed the
  mass the field actually paid, by a bounded and currently unreported amount.
  Quantitative density and carrying-capacity claims therefore remain out of
  scope under `sherwood` — the audit's no-go on that point is not retired by
  choosing it.
- The cost of funded accounting has never been measured. `delivery` gives up the
  agent-metabolism kernel and the chemistry solve while mechanics and toxin
  kernels stay on device, so "we lose the GPU" is directionally right and
  quantitatively unknown.

**Benchmark before committing further.** The GPU cost/benefit is being measured
rather than argued, per uptake mode and per backend, with cost and precision
reported side by side; a speedup quoted without its precision column is not
reportable. The existing `StepProfile` per-phase accumulators are the
instrumentation. `delivery` + GPU is recorded as *blocked* pending a route from
`docs/CUDA_DELIVERY_PARITY.md`, never approximated from a neighbouring arm.

**Self-oracle tests.** The observation that `tests/test_greens_function.cpp`,
`tests/test_fmm.cpp` and `tests/test_octree.cpp` use `concentration_bounded()`
as their own expected value is correct, but those tests are legitimate for their
actual purpose — FMM and superposition accuracy against a direct sum. They were
left alone. The gap was that no test constrained the boundary condition itself,
which is now covered by independent oracles: a wall-gradient residual bound, a
mode-expansion representation of the same boundary-value problem, and a global
mass-balance check.

## Still open — needs an experiment, not an opinion

- QSSA timescale separation, and the structural non-identifiability argument.
  Both are correctly reasoned in the report; both need measurement.
- Whether a clipped reaction should halt a shipped run by default, and the
  magnitude of clipping under `sherwood`. This is the number that says how far
  `sherwood` really is from funded accounting.
- Causal-validation scope. The report's separation of "implements a falsifiable
  coupled EARI/VADI hypothesis" from "has validated EARI/VADI as the causal
  explanation" is adopted. External parameterization, held-out data, competing
  models, and ablations do not yet exist, so claims stay on the former side.

## Round two — the Robin correction design review

The same external platform reviewed the lumen-boundary design before it was
implemented (`gutibm_robin_correction_design_review.md`, review basis PR #359 at
`49dae68`). This section records that review the same way as the audit above:
claim, independent verdict, and what changed as a result. The implementation is
PR #361; the physics and its accuracy budget live in `docs/ROBIN_LUMEN_BC.md`
and are not repeated here.

Context for a reader arriving cold: the lumen at `z_hi` was sealed, i.e. toxin
reflected off it, while `AdvectionField` carries flow through it. That is not a
defensible outflow model, and it is not a small effect — at shipped parameters a
transparent lumen gives 0.35–0.67x and a perfect sink 0.01–0.52x the sealed
toxin concentration, which dwarfs the 9–19% bias from the image-series defect
that PR #359 fixed. The decision taken was the middle treatment, a Robin
mass-transfer condition at the lumen with an impermeable epithelium, because it
is the only one of the three that follows from a physical picture (a well-mixed
lumen bulk above a diffusive boundary layer) rather than from convenience.

### Standing of the review

Substantially correct and materially useful: it changed the implementation in
three places, one of which was a sign error in our own specification that would
have produced a plausible-looking but wrong field. It was also wrong in one
recommendation and unnecessary in two, in each case demonstrably so by
measurement rather than by argument. Its central procedural instruction —
validate the reconstructed total concentration, not the correction — is the most
valuable sentence in it, and it is what exposed the test-design error described
below.

| # | Review point | Verdict |
|---|---|---|
| 1 | Keep `G_Robin = G_sealed + ΔG` and precompute the correction | Adopted, and independently required: a direct mode kernel needs 393 modes at `ρ=1 µm` and 197 at 2 µm, exactly the distances agents evaluate |
| 2 | Expose transfer length `δ`, not `k_c` | Adopted (`toxin.lumen_transfer_length`); `k_c` is not exposed |
| 3 | Do not key tables by receptor field | Adopted. Confirmed independently: ColE1 and ColE2 both target BtuB with different `D_free` and retardation, so one BtuB table would apply the wrong kernel |
| 4 | Declare the state variable and the boundary diffusivity convention | Adopted as a config choice rather than a single declaration — `toxin.lumen_transfer_basis` selects `k_c = D_eff/δ` (default) or `D_free/δ` — because the review's own table is the argument against hard-coding the free basis: shipped pI-driven retardation spans ~55x, so that basis silently turns the lethal-core toxins into an absorbing lumen |
| 5 | Carry normal velocity in the transformed Robin coefficient, and re-derive the sign | **Adopted, and it caught a real error of ours.** See below |
| 6 | Guard near-wall/near-axis with a direct branch | Adopted as a bounded fallback, but the premise was not reproduced — see below |
| 7 | Low screening needs a direct path or a documented rejection gate | Partially adopted: measured, documented, not gated — see below |
| 8 | Advection varies with depth and time, so a table built once goes stale | Confirmed and fixed, though not as ranked — see below |
| 9 | Quantized, bounded cache rather than one table per floating-point tuple | Adopted, and it was load-bearing: mutation moves pI, `retardation_from_pI` moves retardation, so the naive key was unbounded at ~2M Bessel evaluations per miss. Tables are keyed on `(Bi, κH, aH, cutoff/H)` in 1% bins, capped at 64 with LRU eviction and counters in provenance |
| 10 | Validate the reconstructed total, not only `ΔC`; keep positivity a diagnostic, not a clip | Adopted, and it exposed a test-design error of ours — see below |
| 11 | `ρ²` near-axis and `log ρ` beyond a matched radius; `N_s,N_t = 48–64`, `N_ρ = 128–256` | Rejected on measurement — see below |
| 12 | `δ = 100 µm` is an exploratory prior, not a colon default; sweep it; report realized Bi | Adopted, and strengthened: the boundary ships **disabled**, so enabling it is an explicit act, and realized Bi is recorded per cached table in provenance |
| 13 | Record grid, key, tolerance, method, basis and a hash in provenance | Adopted; the hash covers the table values, not the metadata string that describes them |

### Where the review changed the code

**The two walls do not share a boundary law (point 5).** The review's transformed
condition `-D ∂ψ/∂n = (k_c + U_n/2) ψ` is right for the lumen, and following its
instruction to re-derive the sign showed that our own specification had the
epithelium wrong. The epithelium is impermeable, so nothing crosses it and the
correct statement is zero *total* flux, `-D ∂C/∂z + U_z C = 0`; the lumen carries
free advective outflow plus film transfer, `-D ∂C/∂z = k_c C`. In the gauge
variable those give `c_lo = -a` and `c_hi = k_c/D_eff + a` with `a = U_z/(2 D_eff)`
— the drift term enters with *opposite* signs at the two walls. Our spec had used
`+a` at both. The corrected eigenproblem is verified against an independent
reference implementation to 1.6e-10..5.9e-9 across five anchors, and the residual
of each wall law is asserted separately in CI.

**Validating the total, not the correction (point 10).** Following this exposed
that the obvious test is the wrong test. A wall-flux residual computed from the
*interpolated* field is not a check of the boundary condition: a 1e-9 m finite
difference taken inside one 3.1 µm trilinear cell returns that cell's average
slope, so curvature that is invisible in field values (1e-3) reappears amplified
to ~1e-1 in a gradient metric. The boundary laws are therefore asserted on the
direct-mode path (measured 1.8e-4 epithelial, 2.7e-4 lumen), reconstruction is
asserted on field *values* (<= 8.4e-3), and the interpolated wall residual is
retained only as a regression guard whose comment records the decomposition and
states that tightening it requires a finer near-wall grid rather than a looser
tolerance.

**Time-varying flow (point 8).** Confirmed, but the fix is not one of the three
options the review ranked. The staleness that matters is not interpolation error,
it is cache identity: sampling the instantaneous velocity makes the key drift
every step under peristalsis, so every step builds a table and then thrashes the
cache. Tables are therefore built from the peristalsis-free mean profile, with
the full time-varying flow retained in the base image series and in the runtime
drift factor. What that neglects was measured rather than assumed: across
peristaltic factors 0.5–1.5 the omitted change in the correction is 6.0e-3 of the
total field.

### Where the review was wrong, and how we know

**Table coordinates (point 11).** The recommended `ρ²`/`log ρ` scheme and the
48–64 x 128–256 grid were not adopted, because a uniform 33-node `ρ` grid over
the 200 µm cutoff was measured at 1.8e-4..1.8e-3 error against the total field —
already below the dominant error term below. The recommendation would have cost
4–8 MiB per tuple instead of 287 KiB to improve a term that is not the binding
constraint. Grid changes should follow the error map, which is exactly the
review's own instruction.

**Near-axis singularity at the Robin wall (point 6).** The premise was not
reproduced. `|ΔC|/C_Robin` falls to 0.3–1% as `ρ → 0` with both points near the
wall; no logarithmic blow-up appeared. The direct branch is retained as a bounded
safety net, but it is deliberately narrow — a 2048-mode root solve in the hot
loop is a runtime cliff, so widening the predicate to catch "either point near a
wall" was reverted, and direct evaluations are now counted in provenance so any
future cliff is visible in output rather than only in wall time.

**Low screening (point 7).** Real, and milder than implied. Measured
`|ΔC|/C_Robin` rises to 2.4 at `κH = 5e-3`, i.e. roughly 0.4 digits of
cancellation, not catastrophic loss. It is documented as a soft spot with the
direct path reachable, and no accuracy gate was added, because a gate calibrated
on an unmeasured fear would reject usable configurations.

**A subtraction we tried and abandoned.** Defining the correction against the
*image series* rather than against the sealed mode sum would have cancelled the
image base's error, which is attractive since that base is now the dominant term.
It fails: the difference then inherits a `1/ρ`-like singular component and a
uniform `ρ` grid gives 14%–528% error. The sealed mode sum stays as the table's
reference, and the image base's error stays as a documented floor.

### What this round found that neither side had

**The image series does not solve the drift PDE when there is wall-normal flow.**
Translations plus z-reversed reflections satisfy zero diffusive flux at both
walls, but a reflected image with reversed `U_z` solves the mirrored-drift
equation, so the sum solves neither: interior residual 2.9e-1 against 1.4e-3 for
a converged mode sum. It is exact only for wall-parallel flow. At shipped flow
the resulting field error is 1.2e-3..6.9e-3 — sub-percent, and now the *dominant*
term in the Robin accuracy budget, ahead of interpolation. This is inherited from
PR #359, not introduced by the Robin work, and it is recorded as a landmine
rather than fixed here.

**Wall-normal drift in the image series (measured, gated).** Confirmed and quantified. A reflected image with reversed `U_z` solves the operator built from `U'=(U_x,U_y,-U_z)`, leaving exact interior residual `-2 U_z ∂f/∂z`; equivalently the reflected pair enforces zero *diffusive* flux `∂C/∂z=0` rather than the physical sealed law `-D ∂C/∂z + U_z C = 0`, so the sealed slab leaks advective flux `U_z C` at both walls. Measured against a converged sealed eigenmode reference (cross-checked against an independent Hankel closed form to ~1e-10), the relative field error is exactly first order in `Pe_z = U_z H / D_eff`: median `0.114*Pe_z`, worst case `0.439*Pe_z` for `Pe_z <= 0.05`, and machine-precision exact at `Pe_z = 0`. The report's residual metric is not reproducible as a statement about the image sum: a second-order central-difference interior residual normalized by `λC` scores 4.73e1 for the *exact* free-space kernel (true residual zero) versus 4.72e1 for the image sum at `H/500`, and the two stay within 0.2% at `H/2000` and `H/8000`, so that metric is entirely discretization error near the source singularity and cannot distinguish an exact solution from the image sum. Drift-corrected images were derived and rejected on evidence: the exact gauge-variable reflection coefficient `R(γ)=(γ-a)/(γ+a)`, `γ=sqrt(q²+κ²)`, depends on the transverse wavenumber, whereas every per-image exponential weight is wavenumber independent; scanning all constant reflection weights buys at most 1.4-1.5x worst case and does not change the `Pe_z` scaling. The series is therefore gated, not corrected: `qssa.drift_envelope_policy` (warn/error/allow) on a measured envelope `|Pe_z| <= 0.05` (<=2% worst case), with `neumann_drift_envelope_evaluations` provenance. At shipped parameters `profile_alpha=1.5` gives `Pe_z(z_s)=0.0926*(z_s/H)^1.5`; the mid-domain probe is 0.033 (inside the envelope) and the default founder band `z/H in [0,0.5]` carries 6.4e-5..1.6e-3 field error, so no shipped outcome moves and the series arithmetic is unchanged. See `docs/NEUMANN_WALL_NORMAL_DRIFT.md`.

**Enabling the boundary changes biology, not just numbers.** With `δ = 100 µm`,
near-lumen toxin falls to 0.59x (effective basis) or 0.45x (free basis) of the
sealed value, and four toxin tests that assert killing or induction fail. Those
failures are the finding: the lumen treatment, not the image-series fix, decides
whether killing happens at all in those scenarios. The boundary therefore ships
disabled — the default is bit-for-bit the previous sealed result — and flipping
it is a scientific decision that belongs with the GPU/precision benchmark
evidence, not with a library change.
