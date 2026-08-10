# Colony-scale spatial observables

`gut_ibm_tools.colony` clusters agents with a sparse-radius DBSCAN.  The
default radius is twice the median agent-agent nearest-neighbour distance;
this is intended to link typical within-colony neighbours without reaching
the next colony.  `estimate_eps(..., rule="knee")` instead selects the knee of
the sorted k-th-neighbour curve.  Every catalog records the rule diagnostics,
chosen-radius noise fraction, and a 0.5x--4x sensitivity sweep.

The agent table preserves the HDF5 `id` as int64 and includes `colony_id`,
coordinates, `type`, and `n_bi_loci`, so per-kill or per-agent records can be
joined directly on `agent_id`.  The colony table reports size, centroid,
radius of gyration, maximum member radius, nearest colony, producer count and
fraction, and producer-threshold flags for 113, 527, and 1361 producers.

When the genome layer includes `id`, `bi_offset`, and `bi_count`, genotype
means `(type, lineage_id, sorted BI-locus identities)` when lineage is present
and `(type, sorted BI-locus identities)` otherwise.  Each BI-locus identity is
the `(bi_toxin_id, bi_immunity_id)` pair, and the offsets define the slice
`bi_*[bi_offset:bi_offset + bi_count]` for each agent.  This makes purity and
clonality use the recoverable per-agent BI-locus multiset.  Older files lack
these join and offset datasets; they fall back to `(type, lineage_id,
n_bi_loci)` or `(type, n_bi_loci)`, preserving the previous behavior.

## Spatial nulls

`gut_ibm_tools.spatial_stats` returns observed 3-D pair correlation `g(r)`,
Ripley's `K(r)`, and `L(r)-r`, together with both required null ensembles.
The z-stratified positional null fixes total N, labels, per-z-bin counts and
the z marginal, while replacing x/y uniformly over each bin's realized
lateral extent.  It destroys lateral clustering within each z stratum.
The label null fixes positions exactly and permutes labels, preserving every
label count; it destroys spatial segregation associated with those labels.
Both use an injected `numpy.random.Generator`.

Choose the null from the question being asked:

- **Is the point pattern clustered at all?** Call
  `compute_spatial_stats(..., labels=None)` and compare with the
  z-stratified coordinate-randomization envelope.  A whole-box uniform null
  is wrong for this model because the epithelial surface at `z=0` creates
  strong z-structure; it would report expected epithelial layering as
  biological clustering.
- **Are producers more clustered than an arbitrary subset of the same cells?**
  Pass an externally defined producer-status mark (or genotype mark) and use
  the label-permutation envelope.  Positions remain fixed and the marks are
  shuffled, preserving the observed composition.

Never pass spatially derived labels to the label null.  DBSCAN `colony_id` is
the most tempting example and is the wrong choice: colony IDs already encode
spatial proximity, so permuting them is circular and guarantees an apparent
departure.  In the short validation measurement this produced departures at
21--23 of 24 radii, whereas externally defined producer status produced
departures at 0--1 radii per snapshot on the same EARI/VADI data.  The module
does not issue an automatic warning because no general test can distinguish a
legitimate spatially structured mark from an incorrectly supplied
spatially-derived mark without false positives.

The reported curves use a same-mark K/g estimator: only pairs sharing the
same label contribute, while intensity and population size use the whole
realized population.  With one label this reduces exactly to ordinary
unmarked K/g.  A label null is only valid when at least two distinct labels
are present; otherwise the result exposes `label_null_valid=False`, returns
NaN label-null curves, and `summarize_excess(null="label")` raises
`ValueError`.

Intensity is estimated from the realized axis-aligned occupied bounding-box
volume.  These simple estimators apply no edge correction and therefore have
boundary bias.  Observed and randomized curves use exactly the same estimator,
so envelope comparisons remain valid; absolute values near the boundary should
not be interpreted as an unbiased infinite-domain estimate.

The first real-output characterization found no interpretable microcolony-
scale `K(r)` signal at 10--50 µm in the short fixtures: `K(r)` is zero or
numerically tiny there because too few close pairs form.  This is an honest
sparse-population baseline, not a statistics defect.  Those measurements used
160-cell, 900-second EARI/VADI fixtures and one short 600-agent local run;
they say nothing about seven-day campaign behavior.

For a single colony, `nn_colony_id` is `-1` and
`nn_colony_distance` is `NaN`.  An all-noise snapshot returns an empty colony
table with the complete column schema and appropriate NumPy dtypes, rather
than omitting the columns.
