# GAP-02: live packet re-entry (v1)

**Decision:** Require an `Active` destination for either Layer-2 slot. For live pool returns, preserve valid local coordinates. For finite `z >= D`, preserve x/y and use the z-center of the receiving domain's outermost grid cell. This is a numerical insertion convention, not a biological attachment depth. Count reseeds only after insertion; if the packet cannot be placed wholly inside the live box, leave it pooled.

## Grid landing regression

The [focused fixture](../tests/test_layer2_reentry.cpp) exports tag 1 through the live boundary hook and routes it back to slot 0 (segment seed 1234, packet seeds 3 and 4). The new landing assertion failed on `632a509`: at the 30 µm face the returned z was 30 µm, not 29 µm. After the fix:

| Export z (µm) | Receiving z grid | Inserted z (µm) | Next live step |
|---:|---|---:|---|
| 30 | 15 × 2 µm | 29 | Resident |
| 30.5 | 5 × 6 µm (z stride 3) | 27 | Resident |

Each return preserved the exported x/y (7.05556, 11 µm) and identity, counted one reseed, and closed stocks. Residence is observed, not required. The halted states, a failed establishment draw and a one-cell grid with no interior center retained the original packet with zero reseeds. Interior returns and tabulated eligibility were unchanged. These results are for the neutral fixture, not the archived runs.

## Single-slot transfer results before the grid rule

[Neutral fixture](../tests/test_layer2_reentry.cpp): one tagged cell, segment seed 1234, packet seeds 3 and 4. With one transit half-life per routing phase, both streams survive distal exit, attempt reattachment, and draw slot 0. Establishment probability is one. All live steps use a single inactive Fix. Read placement immediately after pool routing; read fate after one live step, **before** another routing opportunity.

| Destination / source | Before (`e5c5421`) | After (`78d170c`) | Next live step (both revisions) |
|---|---|---|---|
| Active live / interior (7, 11, 13) µm | 1 reseed, same position | Same | Same tag resident |
| BloomHalted live / interior | 1 reseed, same position | 0 reseeds; original packet pooled | Not stepped |
| SpatialHalted live / interior | 1 reseed, same position | 0 reseeds; original packet pooled | Not stepped |
| Tabulated / same interior | Active accepts; both halted states retain | Unchanged | Not a live slot |
| Active live / actual boundary export (7.05556, 11, 30) µm | 1 reseed, same position | Same | Same tag washed out to pool |

Stocks closed throughout. The eligibility test checks identity, lineage, position, pool ownership and per-cell reseeds; a retained halted packet's RNG is at the post-target draw, not the pre-attempt state.

**Original failure:** In the 40 × 40 × 30 µm box, the boundary hook exported tag 1 at exactly `z = D`. The accepted cell already failed `0 <= z < D` **before** its subsequent washout. The original witness reported this predicate without asserting it; the grid regression above now checks the landing before the next live step.

## What main does at `28242fd`

| Pool destination | After a packet survives distal exit and draws this target |
|---|---|
| Tabulated | Only `Active` can establish. On success, move cells to `colonists`; otherwise return the packet to the pool. |
| Live | On establishment, move cells to the embedded simulation and increment `reseeds` by cell count, even for `BloomHalted` or `SpatialHalted`. Patch stepping skips both halted states. This transfer can conserve stocks while violating eligibility. |

See [step order](../src/layer2/segment.cpp#L690-L709), [pool routing](../src/layer2/segment.cpp#L716-L763) and [halt conditions](../src/layer2/segment.cpp#L793-L813). The separate Layer-3 arrival path, [`try_deposit`](../src/layer2/segment.cpp#L777-L789), rejects **all** live destinations; it is not pool re-entry.

**Coordinates:** A [patch](../src/layer2/segment.h#L56-L82) has an ID, not a world-space offset. Every [live box](../src/layer2/segment.cpp#L507-L515) starts at zero with common extents `(L,L,D)`, periodic in x/y only. [Tabulated founders](../src/layer2/segment.cpp#L218-L229) start at the local midpoint; [pre-seeded singles](../src/layer2/segment.cpp#L295-L307) start at zero. Live [washout and fragment hooks](../src/layer2/segment.cpp#L570-L605) copy the departing agent. Live insertion resets owner and cached grid index, **not** `Agent::x`.

An in-bounds local point needs no patch-ID translation. But a [boundary export](../src/core/simulation.cpp#L2052-L2091) carries `z >= D`, unlike an imposed washout that can occur below the surface. If accepted unchanged, the former starts at the exit boundary or beyond. [Grid lookup clamps indices](../src/core/domain.cpp#L61-L68), not the stored position. Immediate re-export is a hypothesis, not an observed outcome.

## What was exercised

- **Existing unit coverage (source inspection):** [live influx](../tests/test_layer2_segment.cpp#L604-L630) seeds 60 pool singles, targets the only (live) patch and requires positive reseeds, more than 20 colonists and ledger closure. It stops at first success; it does not check the arriving tags, coordinates or next-step fate. [Handoff](../tests/test_layer2_segment.cpp#L632-L674) checks distal exit with reattachment disabled. [Bloom halt](../tests/test_layer2_segment.cpp#L322-L351) uses tabulated slots, not a transfer into an already halted live slot.
- **Archived measurements:** The [placement archive](../analysis/spec13_placement/out/) (reported build `778c2af`, seeds 20261027–20261031) has 20 incubator runs, each with one live slot among ten. Across their provenance ledgers: 1,569 washout departures, **zero cells reseeded**, all 20 closed. No successful placement to inspect; zero does not count attempts. Four churn ledgers each report one reseeded cell (crypt, seeds 20261029/31, both lysis arms); that segment-wide count does not identify a live destination or landing position.

## Original single-slot test proposal (now exercised above)

Use a deterministic single-slot transfer test with a tagged packet and a known establishment draw. For either halted state, verify the target was selected, no insert or reseed occurred, and the packet remained pooled. On an active target, test an interior input and a boundary export: require exactly one transfer per cell, preserved identity/lineage, finite destination-local coordinates (`0 <= x,y < L; 0 <= z < D`) before the next step, and ledger closure. Record destination status, source and inserted positions, then the tagged cell's next-step fate. An invalid coordinate **on an accepted cell** demonstrates a placement failure; later washout alone does not. Choose an entry rule for invalid points, not a new landing distribution for valid ones.

**Defer multi-slot and restart design.** The [generator](../analysis/spec13_placement/make_configs.py#L76-L92) uses the [scalar parser key](../src/io/input_parser.cpp#L1869) and [one slot initializer](../src/layer2/segment.cpp#L548-L556). The [ten-live reference](../analysis/spec13_placement/REPORT.md#L132-L145) is a proposed study, not an in-tree caller of `audit_patch_indices` or `n_live_patches`; none was found. Multiple slots would also need unambiguous [per-slot provenance](../src/layer2/segment_io.cpp#L230-L233). [Checkpoint write](../src/layer2/segment_io.cpp#L277-L282), [restore](../src/layer2/segment_io.cpp#L594-L598) and [Layer 3](../src/layer3/chain.cpp#L99-L103) explicitly refuse live slots. Keep those refusals; no checkpoint redesign, all-live comparison or population outcome is needed to decide the re-entry contract.
