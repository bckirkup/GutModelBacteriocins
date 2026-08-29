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
| 3 | Configuration fails open (missing file, malformed JSON, unknown keys) | Confirmed. A whitelist of `ConfigError` messages does rethrow, so it is not uniformly fail-open, but the paths that matter most are |
| 4 | GPU per-species diffusion completion suppresses CPU fallback | Confirmed. `apply_diffusion` returns "any species succeeded" and the pipeline falls back only when that aggregate is false. Their specific `nz=1025` reproducer was not checked |
| 5 | GPU mechanics omits CDI corpses | Confirmed **and understated**. `cdi.enabled` defaults true with 300 s persistence, so this is not a narrow guard: every GPU run in which any cell dies diverges from CPU for five simulated minutes per death |
| 6 | Checkpoint restart does not restore RNG state | Confirmed, low novelty — `docs/BRANCHING_FROM_CHECKPOINTS.md` already states a fork is a population-state continuation, not a bit-identical restart |
| 7 | Post-division SOS hazard applies to both parent and daughter | Confirmed: 1.98% per division event against a documented 1% "per division". An estimand ambiguity, cheap to settle |
| 8 | Comet-tail analysis fabricates geometry | Confirmed, and the most damaging analysis defect. Coordinates are assigned by `linspace` over flattened storage order, so every comet-tail number produced to date is a function of HDF5 layout, not of flow |
| 9 | Resident-retention estimator measures lineage-label persistence | Confirmed |
| 10 | Flagship diversity scenario has no immigration block and no grid output | Confirmed, both halves, against a README that describes periodic lumen immigration |
| 11 | Taylor–Aris dispersion toggle is unused | Confirmed: zero callers, defaulted on, documented as active |
| 12 | VBF drag and carrying capacity unused | Half wrong, right conclusion — see below |
| 13 | Requested HDF5 output can fail open | Confirmed: on create/validate failure the file is removed and the run proceeds |
| 14 | Provenance and compiler contract incomplete | Confirmed: git SHA falls back to `unknown-git-unavailable`, and `<format>` is used widely while the README never states the GCC 13+ requirement |
| 15 | Python manifest writes are non-atomic | Confirmed: destination opened in `"w"`, no temp file, no fsync, no retained generation |
| 16 | Ethanolamine absolute units off by 1000 | Confirmed, already recorded in `docs/UNITS_AUDIT.md`. `eut_km` is off by the same factor, so the Monod penalty is numerically unchanged — labels and coupling only |

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
