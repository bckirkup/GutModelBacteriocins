# T2 Colony-Challenge Measurement: finite-colony toxin efficacy

## Executive finding

The controlled T2 assay measures substantially less target occupancy than the
earlier point-source arithmetic predicted.  The measured half-occupancy
thresholds are:

| Target distance from nearest producer | Measured threshold |
|---|---:|
| 10 µm | 300–1000 producers |
| 30 µm | 1000–3000 producers |
| 50 µm | above 3000 producers |

These are higher than the earlier analytic predictions of approximately 113,
527, and 1361 producers at 10, 30, and 50 µm respectively.  The discrepancy
is geometric rather than a failure of the toxin solver: the earlier
calculation treated the colony as a point source.  A real colony has physical
extent, so most producers are farther from the target than its nearest cell.
That correction grows with producer count, exactly where the threshold
matters.

This document is the result of a controlled mechanism assay.  It does not
establish campaign-level killing or retention.

## Experimental design

The assay used a non-periodic `220 × 150 × 150 µm` domain with a 5 µm chemical
grid, 60 s biology steps, and a 600 s total run.  Advection, motility,
mechanics, CDI, and crypt migration were disabled.  Each case contained an
induced producer colony and one target placed at a requested distance from the
nearest producer.

Producer counts were:

```text
N = 1, 10, 30, 100, 300, 1000, 3000
```

Target distances were:

```text
d = 5, 10, 20, 30, 50, 100 µm
```

Producer centers were placed on a centered cubic lattice with 1.2 µm
cell-center spacing.  The first `N` lattice points were selected in increasing
distance from the colony center.  The resulting approximate outer radii
were:

| N | Outer radius |
|---:|---:|
| 1 | 0.00 µm |
| 10 | 1.70 µm |
| 30 | 2.40 µm |
| 100 | 3.60 µm |
| 300 | 4.95 µm |
| 1000 | 7.40 µm |
| 3000 | 10.80 µm |

Distance is always measured from the nearest producer, not from the colony
centroid.  Consequently, at `N=3000` the colony radius is comparable with the
shortest target distances; treating this colony as a point source is not a
small approximation.

The assay used the model-default colicin kill rate:

```text
kill_rate_colicin = 1e-3 s^-1
```

The target was run as a non-lethal observer.  It remained present for the
complete exposure while the toxin field was sampled at its assigned grid
cell.  B12 competition was included using the model's apparent Kd:

```text
apparent_Kd =
    5e-7 × (1 + 1e-3 / 1e-6)
  = 5.005e-4 mol/m³
```

For each post-lysis step, competitive occupancy and hazard were evaluated as:

```text
occupancy = C / (5.005e-4 + C)
hazard = 1e-3 × occupancy × 60
```

The measured probability was then calculated without a Bernoulli draw:

```text
p_kill = 1 - exp(-Σ hazard_t)
```

This is the transferable method.  The field is deterministic and the
per-step hazard is already known, so sampling a stochastic death process is
unnecessary when the quantity of interest is its analytic probability.

The nominal design was therefore 42 observer runs rather than 840 uniformly
replicated Monte Carlo runs.  It ran on two local CPU cores without a cluster.
Three nominal cells (`N >= 300, d=100 µm`) had no valid target because
immigration failed; they are reported as invalid, never as biological zeros.

The sweep harness and raw HDF5 artifacts remain untracked under `t2_assay/`;
they are not part of this documentation change.

## Deterministic surface

The table reports:

```text
analytic p_kill [mean competitive occupancy over the nine exposure steps]
```

| N | 5 µm | 10 µm | 20 µm | 30 µm | 50 µm | 100 µm |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.005 [0.010] | 0.002 [0.005] | 0.001 [0.002] | 0.0005 [0.001] | 0.0002 [0.0004] | 0.00005 [0.00009] |
| 10 | 0.044 [0.083] | 0.017 [0.032] | 0.007 [0.013] | 0.004 [0.008] | 0.002 [0.003] | 0.0003 [0.0006] |
| 30 | 0.062 [0.119] | 0.039 [0.073] | 0.021 [0.038] | 0.013 [0.024] | 0.005 [0.009] | 0.001 [0.002] |
| 100 | 0.150 [0.301] | 0.104 [0.203] | 0.057 [0.108] | 0.034 [0.063] | 0.015 [0.028] | 0.004 [0.008] |
| 300 | 0.224 [0.470] | 0.173 [0.352] | 0.110 [0.217] | 0.087 [0.169] | 0.039 [0.073] | invalid |
| 1000 | 0.299 [0.659] | 0.299 [0.659] | 0.228 [0.478] | 0.170 [0.345] | 0.092 [0.180] | invalid |
| 3000 | 0.367 [0.846] | 0.340 [0.768] | 0.311 [0.689] | 0.252 [0.537] | 0.184 [0.377] | invalid |

The `N=1000` values at 5 and 10 µm are identical because both target
positions map to the same 5 µm chemical grid cell.  This is a resolution
limit, not evidence that the underlying continuous field is identical at
those two distances.

The occupancy transitions in the table give the measured thresholds:

- at 10 µm, occupancy crosses 0.5 between `N=300` and `N=1000`;
- at 30 µm, it crosses 0.5 between `N=1000` and `N=3000`;
- at 50 µm, it remains below 0.5 at `N=3000`.

The point-source predictions of 113, 527, and 1361 producers therefore
underestimate the finite-colony producer counts needed for half occupancy.

## Why the measurement is credible

Three independent checks support the interpretation.

### Direct QSSA and FMM agree

At `N=300`, `d=10 µm`, the direct-QSSA and FMM toxin grids agreed at all
11 saved snapshots:

```text
maximum absolute difference = 0
maximum relative difference = 0
```

The FMM path also completed the `N=3000` cases at distances 5–50 µm after the
degenerate-cluster guard was merged.

### Observer hazard matches killed-target provenance

An independent killed-target run at `N=300`, `d=5 µm` recorded a first
provenance hazard of:

```text
0.03969635391883309
```

The observer run reconstructed the first hazard from the target-cell toxin
trace as:

```text
0.039696353918833166
```

The difference is floating-point roundoff.  This verifies that the observer
trace is using the same competitive occupancy and hazard calculation as the
kill path.

### Monte Carlo comparison

Existing susceptible-target runs were compared with the observer-derived
probabilities.  The earlier 20-replicate `N=100` rows are consistent at every
distance.  The exact one-sided binomial tail probabilities are shown below;
they quantify how surprising the observed count is under the deterministic
probability.

| N | d | Monte Carlo kills | Trace p_kill | One-sided tail |
|---:|---:|---:|---:|---:|
| 100 | 5 µm | 2/20 | 0.150 | 0.405 |
| 100 | 10 µm | 1/20 | 0.104 | 0.370 |
| 100 | 20 µm | 0/20 | 0.057 | 0.311 |
| 100 | 30 µm | 0/20 | 0.034 | 0.504 |
| 100 | 50 µm | 0/20 | 0.015 | 0.741 |
| 100 | 100 µm | 0/20 | 0.004 | 0.918 |
| 300 | 5 µm | 4/5 | 0.224 | **0.0103** |
| 300 | 10 µm | 2/5 | 0.173 | 0.209 |
| 300 | 20 µm | 1/5 | 0.110 | 0.443 |
| 300 | 30 µm | 1/5 | 0.087 | 0.367 |
| 300 | 50 µm | 1/5 | 0.039 | 0.178 |
| 1000 | 5 µm | 0/3 | 0.299 | 0.344 |
| 1000 | 10 µm | 0/3 | 0.299 | 0.344 |
| 1000 | 20 µm | 0/3 | 0.228 | 0.461 |
| 1000 | 30 µm | 0/3 | 0.170 | 0.572 |
| 1000 | 50 µm | 0/3 | 0.092 | 0.748 |

The `N=300, d=5 µm` result is an unexplained outlier: 4 of 5 targets died
against a trace-derived probability of 0.224.  Its one-sided binomial tail is
0.0103.  It is retained here rather than explained away or used to tune the
analytic result.  The remaining comparisons, especially the 20-replicate
`N=100` rows, support the deterministic collapse but do not replace a larger
Monte Carlo campaign if that outlier becomes scientifically important.

## Limits and interpretation

### This is a mechanism assay, not a campaign

The controlled synthetic assay characterises the toxin mechanism under fixed
geometry.  It does not establish that a host-scale campaign reaches these
producer counts, maintains these colonies, or produces these killing
probabilities.

The assay delivers nine post-lysis exposure steps.  Its per-run ceiling is
therefore:

```text
P_max = 1 - exp(-9 × 1e-3 × 60) = 0.417
```

The earlier `0.259` value came from the mistaken five-step assumption:

```text
1 - exp(-5 × 1e-3 × 60) = 0.259
```

It is not the ceiling of this harness and must not be interpreted as a
biological plateau.

### Grid resolution limits short-range interpretation

The 5 µm grid aliases the `N=1000` cases at 5 and 10 µm.  Target placement is
continuous, but the recorded toxin concentration and reconstructed occupancy
are read from the assigned grid cell.  The short-range surface therefore has
finite spatial resolution and should not be interpreted as a continuous
sub-grid measurement.

### Realized spatial structure remains the next measurement

The available spatial-statistics runs have `K(r)` identically zero at
10–50 µm.  In those runs there are no close pairs at those radii.  This means
the measured thresholds are demanding in a system whose realized structure is
not obviously producing large, tight, close colonies.

It does **not** show that campaigns fail to reach those structures; that is
untested.  It makes the question measurable: future burn-in and fork
experiments should record the colony catalog and compare realized
`K(r)`, colony sizes, and producer counts against this finite-colony
threshold surface.  The spatial-observable conventions and null selection are
documented in [`COLONY_OBSERVABLES.md`](COLONY_OBSERVABLES.md).

### Immigration validity is part of the measurement

The `N=3000, d=100 µm` placement failed consistently in the original assay.
That was a search deficiency in the immigration proposal, not a biological
zero, and it is now handled by the separate immigration fix.  Invalid
placement must continue to be reported explicitly in future sweeps.

## Relation to the experimental ladder

This assay is an L1 microcolony mechanism measurement in the terminology of
[`MULTI_SCALE_EXPERIMENTATION.md`](MULTI_SCALE_EXPERIMENTATION.md).  It should
feed later burn-in/fork experiments rather than be promoted directly into a
host-scale claim:

1. measure realized colony sizes and `K(r)` at the relevant radii;
2. identify whether the campaign produces colonies in the 300–3000 producer
   range;
3. fork those realized states into controlled toxin challenges;
4. reserve campaign-scale runs for retention and sweep questions after the
   local mechanism is characterized.
