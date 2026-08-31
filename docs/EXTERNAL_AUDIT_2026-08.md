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
| 1 | Two-wall Neumann image series incomplete and duplicated | Confirmed, and worse than described — see below |
| 2 | `uptake_limit` ships as `none`, so growth is unfunded by realized removal | Confirmed factually; the remedy is a decision, not a patch — see below |
| 3 | Configuration fails open (missing file, malformed JSON, unknown keys) | **Confirmed and fixed.** JSON-vs-legacy selection now uses file shape rather than an error-message whitelist: once a root JSON object is found, every parse failure is fatal; files without one retain legacy flat-key parsing. Named-file open failures, malformed arrays, unknown keys, and invalid known values fail closed by default. `GUTIBM_STRICT_CONFIG` is inverted: unset/`1` is strict, while `0` is the explicit lenient escape hatch. |
| 4 | GPU per-species diffusion completion suppresses CPU fallback | **Confirmed and fixed.** A whole-field pre-flight now declines GPU diffusion before any launch when an eligible species exceeds the line cap, so the pipeline's existing host fallback covers every species. The mixed-eligibility window is exactly `nz == 1025`: Dirichlet solves `nz - 1`, Robin/Flux solve `nz`; above the cap every species is refused uniformly. Their `nz = 1025` reproducer is covered host-side for the predicates and on the T4 for the device path; a species the host does not diffuse never forces a fallback |
| 5 | GPU mechanics omits CDI corpses | **Confirmed and fixed.** Host and device now share the same death-time/persistence predicate; fresh CDI corpses remain mechanical obstacles while expired corpses are excluded. Host-side coverage is local, while the device assertion is evidenced only by the T4 CI job |
| 6 | Checkpoint restart does not restore RNG state | Confirmed, low novelty — `docs/BRANCHING_FROM_CHECKPOINTS.md` already states a fork is a population-state continuation, not a bit-identical restart |
| 7 | Post-division SOS hazard applies to both parent and daughter | Confirmed: 1.98% per division event against a documented 1% "per division". An estimand ambiguity, cheap to settle |
| 8 | Comet-tail analysis fabricates geometry | **Confirmed and fixed.** Validation now uses physical HDF5 cell centers in `(z, y, x)` flattening order, the model's positive-x distal-flow convention, and an explicit producer reference with periodic minimum-image wrapping |
| 9 | Resident-retention estimator measures lineage-label persistence | Confirmed |
| 10 | Flagship diversity scenario has no immigration block and no grid output | **Confirmed and fixed, with a substantive negative result.** Both halves stood: no `immigration` block and `grid: 0`. The example now ships continuous luminal immigration and 6-hourly named-species grid output, and an out-of-domain `z_slab` band now fails closed instead of silently clamping immigrants onto a domain face. Enabling the documented mechanism does **not** change the flagship regime at shipped parameters — see below. |
| 11 | Taylor–Aris dispersion toggle is unused | **Confirmed and resolved by removal.** The dead, default-on toggle and uncalled implementation were removed rather than wired into the isotropic transport model. At `z=h`, the measured `D_taylor/D` enhancements are `6.379e-9` for oxygen (`D=2e-9`), `6.379e-5` for carbon (`D=2e-11`), `4.023e-2` for shipped ColE1-like colicin (`D=7.964e-13`), and `2.552` for a hypothetical `D=1e-13`; the enhancement falls as `(z/h)^3`. The long-time limit is valid for the shipped colicin (`t/t_diff=3.44`) but not for the large-effect hypothetical case (`0.43`). The `210` prefactor is also the Poiseuille result while the shipped profile exponent is `1.5`. Honest wiring would require an anisotropic streamwise transport kernel and re-derivation of both the image series and Robin correction table. Reinstating it requires that work. |
| 12 | VBF drag and carrying capacity unused | Half wrong, right conclusion; **resolved by removal** — see below |
| 13 | Requested HDF5 output can fail open | **Confirmed and fixed.** Requested HDF5 output now throws before compute when path validation or file creation fails; under MPI, rank 0 makes the decision and broadcasts it so every rank fails together rather than continuing without the requested record. The failed output file is no longer silently removed, and the error names the file plus the underlying reason where available. |
| 14 | Provenance and compiler contract incomplete | Confirmed: git SHA falls back to `unknown-git-unavailable`, and `<format>` is used widely while the README never states the GCC 13+ requirement |
| 15 | Python manifest writes are non-atomic | **Confirmed and fixed.** Both run-level JSON writers (`batch_manifest.save_manifest`, `batch_runner._write_json_output`) now route through the hardened `path_utils` layer, which serializes fully before opening anything, writes a same-directory temp file, `flush()` + `os.fsync()`, `os.replace()`, then a best-effort directory fsync; the temp file is removed on failure. Regression tests interrupt the serializer and assert the previous generation is byte-identical, still parseable, and that no temp file is left behind — they fail on the pre-change code with an empty destination (`assert b'' == b'{\n  "versi...'`). |
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

### Claim 10 measurements and resolution

The flagship geometry is not runnable to its 7-day horizon here, so both arms
were measured on a scaled fixture: 400×400×100 µm (x and y reduced 5×), all
biology, `grid_dx`, mucus thickness, `radial_turnover` and `distal_transit`
unchanged, founders scaled to the same areal density as the shipped 500 + 100
in 2×2 mm (20 residents + 4 immigrants), seed 42, 12 h = 720 steps of 60 s,
serial build, `termination_cause_code = horizon_reached`. The two arms differ
only in `immigration.enabled`.

| Observable at horizon | immigration on | immigration off |
|---|---:|---:|
| Immigration events (cumulative) | 14 | 0 |
| Live immigrants (type 2) | 0 | 0 |
| Live residents (type 1) | 121 | 195 |
| Boundary outflow (cumulative) | 61 | 33 |
| SOS lysis deaths (cumulative) | 4 | 3 |
| `monochromatic_score` | 1.000 | 1.000 |
| Grid dumps emitted (3 species each) | 13 | 13 |

1. **The mechanism is reachable and fires.** 14 Poisson events at the shipped
   1/h rate, 11 distinct post-founder immigrant tags observed in the agent
   snapshots, first-seen z between 84.7 and 99.4 µm.
2. **Immigrants do not establish.** Observed residence was 0–4800 s (median
   600 s, the snapshot interval). Radial advection reaches the lumen boundary
   from the arrival band in a few hundred seconds at `profile_alpha = 1.5`
   and `radial_turnover = 5400 s`, against a resident doubling time of
   `ln2/5.0e-4 ≈ 1390 s`. The washout interaction is real and was measured,
   not assumed: cumulative boundary outflow nearly doubled (61 vs 33) while
   `outflow_washout` stayed 0 in both arms (`washout.trap` is `emergent`).
3. **So the headline observable does not discriminate.** Type-based resident
   retention is 100% in both arms — the initial immigrant cohort dies out
   without any immigration — and `monochromatic_score` is degenerate at 1.000
   once one type remains. Enabling the documented mechanism changes the influx,
   not the outcome. What it buys is that the washout half of the Advective
   Double-Bind is now exercised and recorded instead of being claimed only in
   prose. Whether immigrants *should* establish at these parameters is a
   scientific question about `mu_max`, the luminal arrival depth and the
   turnover time; no parameter was retuned to make the example look better.
4. **Grid output exposes a second limitation.** With grid dumps finally
   emitted, the bacteriocin fields are identically zero in most of them:
   release is tied to SOS lysis (3–4 events per 12 h arm) and the field decays
   within a few steps, so a strided grid dump is unlikely to catch an
   induction. The comet-tail and lethal-halo claims need a dedicated
   high-cadence grid window, which is why the shipped stride is 6 h rather
   than hourly — three named species at the flagship grid size are already
   ~5 GB per week of simulated time.
5. **MPI equivalence.** Immigration decisions are drawn from the replicated
   stream (`seed ^ 0x9e3779b97f4a7c15`) and are rank-count invariant: a 2 h
   control fired exactly 1 event at both 1 and 2 ranks, with globally unique
   agent tags and no duplicated insertion (only the owning rank constructs the
   cell). Serial and 2-rank *trajectories* do diverge, but they diverge
   identically with immigration disabled (5 vs 11 live agents in the control),
   so that divergence is the pre-existing rank-partitioned agent RNG, not
   something immigration introduced. RNG consumption is unchanged when
   immigration is disabled, and the claim-6 checkpoint RNG-state limitation
   is not worsened.
6. **Not verified locally.** The full 7-day 2×2 mm horizon (grid alone is
   ~5×10⁷ cells) and any GPU path.

## Where we corrected the report

**#12 VBF drag — right conclusion, wrong mechanism.** `drag_force()` *is* called
from the physics module, so "stored but unused" is inaccurate. The real defect is
different and slightly worse: `a.v` is written only by that drag term, drag is
`-c·v`, and mechanics is overdamped and never writes `a.v`, so a velocity that
starts at zero stays identically zero and the `a.x += a.v·dt` update is a
permanent no-op. `local_capacity()` genuinely has no callers, so the
carrying-capacity half stands as written.

Measured, not argued: over a 60-step VBF-enabled serial run (60 µm × 60 µm ×
40 µm domain, 200 founders, `bio_dt = 60` s), `max|a.v|` over every agent and
every step is exactly `0`, and so is `max|a.v|·dt`. The instrumentation read
`a.v` off the agents directly rather than any flux counter, so it is not
exposed to the `*_step`-versus-`*_last_step` accounting trap.

The channel was also unusable as written, which decided the fix. With
`drag_coeff = 1e-9` N·s/m, a measured agent mass of `9.13e-16` kg and
`dt = 60` s, the explicit update `v ← v(1 - c·dt/m)` has `c·dt/m = 6.6e7`: had
`v` ever been nonzero it would have diverged immediately. Making drag entrain
agents toward the local fluid velocity — the physically defensible alternative —
therefore needs an implicit or analytic relaxation at `bio_dt`, not the existing
explicit step, and it would change trajectories and so every spatial observable.
We took removal instead: the velocity integration, the `a.x += a.v·dt` update,
`drag_force()`, `local_capacity()`, and the dead `vbf.drag_coeff` /
`vbf.carrying_cap` config keys are gone, mirroring the Taylor–Aris decision
(claim 11) to delete a documented-but-dead path rather than wire it silently.
The change is numerically inert: the removed update added exactly `0` to every
position. `Agent::v` is kept, because MPI transfer and checkpoint layouts
serialize it, and it remains zero. Mechanics is now documented as purely
overdamped with translation from advection plus mechanics displacement.

The accompanying test asserts the invariant (`a.v` exactly zero for every agent
at every step) and is a characterization test: it passes both before and after,
which is the point — the removal changed no number. The fails-before-passes-after
tests in this PR are the Python atomicity ones (claim 15) and the mucin
dimensional one.

**Daughter placement near bounded walls (Section 4).** `Domain::apply_pbc`
clamps non-periodic axes, so near-wall daughters are never lost from the domain —
they accumulate on the wall. The artifact is a surface density bias, not
immediate loss, which makes the risk milder and the test cheaper than claimed.

## What we found that the report did not

**The dynamic mucin liberation prefactor is dimensionally wrong, and it is the
parameter that is wrong, not the expression.** This was raised as a suspicion
against `dynamic_mucin_liberation()` and checked rather than assumed. The return
value feeds `chem.reac`, so it must be `mol m^-3 s^-1`, and the static
alternative it replaces (`vbf.mucin_liberation = 5e-5 mol m^-3 s^-1`) fixes those
dimensions for both branches. With `vbf.density` in cells m^-3, the prefactor
must be `mol cell^-1 s^-1`, but `k_liberation` was declared `1/s` with default
`1e-4`, giving `1e-4 * 1e11 * 0.909 = 9.1e6` against `5.0e-5`, a factor `1.8e11`.
The expression is correct once `k_liberation` is read as a per-cell specific
rate, which is also the intended physics, so the fix is to redeclare it
`mol cell^-1 s^-1` with default `5.0e-16` — calibrated so the dynamic path
reproduces the static term at the default background density — and to change no
host or device code, which is why the backends cannot diverge here. `mucin.enabled`
is `false` in every shipped example and default, so no published-style number
moves. Full algebra in `docs/UNITS_AUDIT.md` §11; the unresolved half (liberation
is not limited by the mucin present, and ignores the z-gradient weighting) stays
open in the landmine table.

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

**Cross-boundary conjugation over-counted rather than lost.** The standing
landmine row described a cross-rank plasmid transfer as lost when ghosts are
cleared. Measurement says otherwise: `FixConjugation` had no ghost-recipient
guard, so a boundary-straddling pair was attempted **twice** per step — once on
the donor's rank, where the recipient is a ghost, and once on the recipient's
owning rank, where the donor is a ghost. Only the second write survives the
ghost clear, but the first still increments `conjugation_transfers` and still
appends a lineage HGT edge. A deterministic fixture (three pairs straddling the
4-rank slab boundaries at `x = 25/50/75 µm` plus one interior control pair,
`p_transfer` forced to 1) recorded 7 transfers and 7 HGT edges for 4 actual
acquisitions at `nprocs = 4`, against 4/4/4 in serial over identical geometry.
Plasmid spread itself was unbiased: with 600 boundary-straddling and 600
interior pairs at `p = 0.3`, cross-boundary recipients acquired at `p`, not at
`1-(1-p)^2`, because the donor-side duplicate trial is what gets discarded. The
defect was therefore in the event ledger and the recorded genealogy, not in the
biology — every phantom HGT edge named a recipient whose genome never changed.
The fix commits acquisitions only on the recipient's owning rank; ghost donors
remain valid, so each pair is attempted exactly once globally with no new
communication.

**Scale of the affected fraction at `nprocs = 4`.** Sampled from serial runs of
the shipped examples (chemistry grid coarsened to fit a 7.5 GB box; pair
geometry is unaffected), classifying every donor/recipient pair within the 4 µm
pili reach by the cell-aligned 4-slab x decomposition: `eari_vadi_validation`
(`Lx = 500 µm`) 10 of 3368 pairs cross a boundary (0.30%), `single_colony`
(`Lx = 1000 µm`) 2 of 2266 (0.09%), `diversity_paradox` (`Lx = 2000 µm`) 0 of
264. Mean `|Δx|` over in-reach pairs is 1.41 µm, close to the 3r/8 = 1.5 µm
uniform-in-ball value, so `nprocs·E|Δx|/Lx` predicts 0.6% at `Lx = 1000 µm` and
0.3% at `Lx = 2000 µm`; the measured fractions sit at or below that because
microcolonies are clustered and only occasionally sit on a boundary — the
occurrence is bursty, one sample reaching 3.2%. So multi-rank HGT counts at
`nprocs = 4` were over-reported by well under 1%, scaling with rank count and
inversely with slab width. Separately, no conjugation-eligible pair (a
plasmid-free recipient within pili reach of a plasmid-bearing donor) occurred at
all in any sampled window of any shipped example: `single_colony` ships only the
conjugative strain, and in the two-strain examples microcolonies are clonal, so
the plasmid-free strain never came within 4 µm of a donor. Under shipped
defaults the over-count is unreachable, which is the reason no published
multi-rank result is at risk. `domain.ghost_width` (10 µm default) must stay at
or above the pili reach (4 µm default) for the recipient's owner to see the
donor at all; `FixConjugation::init` now warns when it does not, since that is
the one configuration where a genuine loss remains possible.

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

**Prescribed oxygen mass was not removed in every directional solve.** In the
default Dirichlet delivery mode, the z solve omitted prescribed mass, so the
old `fractions[0] - fractions[3] > 0.5` expectation was reachable only because
respiration was funded from oxygen that the field had not removed. At the
highest tested rung, the pre-fix probe credited `3.958e-18` mol while the
field-side removal was `1.002e-18` mol; after the prescribed-mass fix, those
values were `1.227e-19` mol credited versus `9.927e-19` mol removed. The
post-fix funded fractions for oxygen concentrations
`[0, 1.0e-6, 2.5e-6, 5.0e-6]` mol/m^3 are
`[1.000000000, 0.993979919, 0.984949798, 0.969899597]`.

The user approved replacing the obsolete magnitude expectation with the
measured band (`fractions[0] - fractions[3] > 0.02` and
`fractions[3] < 0.98`) and with the invariant that credited respiration over a
step does not exceed the oxygen inventory decrease plus boundary influx.

Measured at the shipped epithelial oxygen default of `55e-6`, funded
respiration is delivery-limited to approximately one third of demand at
`2e-6` voxel resolution. Saturation is reached near `166e-6` and is clean at
`200e-6`; this is a real constraint on the funded-uptake default decision, not
a test artifact. The linear solve is linear in the prescribed mass, so the
delivery limit lifts proportionally with oxygen availability.

**The reported unaccounted oxygen removal was a false alarm.** The per-step
oxygen ledger closes to machine precision when normalized by gross channel
traffic: the measured worst signed residual ratios were `-3.5e-13` (uniform
initialization, gradient on), `1.6e-13` (model initialization, gradient on),
`-3.0e-12` (uniform initialization, gradient off), `-3.0e-12` (model
initialization, gradient off), `1.6e-13` (model initialization, gradient on,
10 steps), and `-5.7e-12` (uniform initialization, gradient off, 10 steps).
The symptom came from reading `boundary_step` after
`commit_boundary_and_reaction_step()` had zeroed it; the correct
`boundary_last_step` value for the reported arm was `-9.57e-18` mol. The
probe also overwrote the field with a uniform `55e-6` mol/m³, far above the
model's depth-decaying reference profile, producing a `-49.6%` first-step
inventory change; using the model initialization instead produced `+0.53%`.
There is no missing oxygen removal channel.

The same measured probe found that the shipped-default delivery limit is
driven by the z-gradient, not voxel size alone. At the agent depth
(`iz = 7`), local oxygen is `3.02e-5` mol/m³ versus the `5.5e-5` epithelial
value. With the gradient enabled, rationing fires, delivery reduction is
`2.75e-18` mol, and funded respiration is approximately `0.33` of demand;
with the gradient disabled, rationing never fires and full demand
(`4.08e-18` mol) is funded. Raising `epithelial_conc` therefore does not lift
the limit while the gradient remains active.

The oxygen ledger has two removal representations depending on whether delivery
is enabled. My first closure identity double-subtracted the delivery VBF mass:
`vbf_sink` is a reporting share contained in `total_sink_realized` for delivery
species, while it is an independent channel for non-delivery species. The
seventh regression arm exposed this accounting-representation hazard; at the
shipped `vbf_sink = 1e-3` and delivery enabled, the pre-fix clamped ledger
reported a VBF share of approximately `6.18e-19` mol/step against approximately
`1.95e-18` mol/step total; after signed accounting, the same arm reports
approximately `6.01e-19` mol/step. The seventh arm's closure is
`1.79e-12` gross-relative after the fix, against `5.47e-3` before it. This is
not a mass-conservation defect: the field removes the VBF mass exactly once;
the difference in reported share is evidence that the clamp moved the ledger.

The corrected seventh arm then exposed a separate clamped realized-sink defect.
The A–E matrix exonerated the z-gradient and localized the residual to the
agent-consuming delivery path. Directional instrumentation localized it to the
periodic x and y solves; the z solve was exact. The replicated periodic clamps
dropped `2.62e-20` mol in x and `3.16e-20` mol in y, `5.785e-20` mol total,
matching the `5.47e-3` gross-relative residual to eight digits. At full demand,
2500 owned cells went negative mid-step, with a most-negative excursion of
`-3.37e-3` mol/m³, roughly sixty times the ambient field value. The sink
created `1.086%` of the step's total realized removal in those cells. End-of-step
positivity still passed because the epithelial Dirichlet face refilled them, so
the #338 rationing loop never fired: its positivity guarantee is end-of-step
only and does not constrain intermediate directional states. Whether to enforce
intra-step positivity remains an open decision; this implementation uses signed
ledger accounting and does not settle that question.

**Non-delivery GPU oxygen VBF remains a discretization divergence.** Delivery
enabled oxygen is a first-order sink in the implicit Route B diagonal, while
ordinary non-delivery GPU VBF intentionally retains the shipped explicit
reaction update. For the shipped `oxygen.vbf_sink = 1e-3 1/s` and `dt = 60 s`,
the explicit concentration factor is `1 - s·dt = 0.94`, versus the implicit
factor `1/(1 + s·dt) = 0.9433962264`. The one-step removed fractions therefore
differ by `6.0%` relative to implicit removal; the resulting field differs by
`0.0033962264·C` in absolute concentration, or `0.36%` relative to the
implicit concentration. This is an open host/device discretization divergence
on the non-delivery GPU VBF path. Its scientific significance is left for a
separate decision.

**The image series does not solve the drift PDE when there is wall-normal flow.**
Translations plus z-reversed reflections satisfy zero diffusive flux at both
walls, but a reflected image with reversed `U_z` solves the mirrored-drift
equation, so the sum solves neither: interior residual 2.9e-1 against 1.4e-3 for
a converged mode sum. It is exact only for wall-parallel flow. At shipped flow
the resulting field error is 1.2e-3..6.9e-3 — sub-percent, and now the *dominant*
term in the Robin accuracy budget, ahead of interpolation. This is inherited from
PR #359, not introduced by the Robin work, and it is recorded as a landmine
rather than fixed here.

**Enabling the boundary changes biology, not just numbers.** With `δ = 100 µm`,
near-lumen toxin falls to 0.59x (effective basis) or 0.45x (free basis) of the
sealed value, and four toxin tests that assert killing or induction fail. Those
failures are the finding: the lumen treatment, not the image-series fix, decides
whether killing happens at all in those scenarios. The boundary therefore ships
disabled — the default is bit-for-bit the previous sealed result — and flipping
it is a scientific decision that belongs with the GPU/precision benchmark
evidence, not with a library change.
