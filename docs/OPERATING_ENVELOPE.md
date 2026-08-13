# GutIBM operating envelope

This document defines the model-side operating envelope for discrete
Enterobacteriaceae in the simulated colonic mucus slab. The companion
configuration inventory, `/home/ubuntu/envelope_geometry.md`,
contains the per-config geometry and initial-density calculations.

## 1. Unit convention and conversion

`Simulation::dysbiosis_threshold_exceeded()` compares the global live-agent
count divided by the **global domain volume in mL** against
`dysbiosis_threshold` in cells/mL. The model domain is a mucus slab: x-y is
mucosal footprint and z is slab thickness.

The primary operating-envelope unit is therefore:

```text
rho_vol = global live agents / mucus-slab volume [cells/mL]
```

The literature reports mucosal bacteria as bacteria/mL of mucus, CFU per
biopsy, or CFU/g. These are not interchangeable without geometry and mass
assumptions:

| Reported as | Typical method | Conversion to model units |
|---|---|---|
| bacteria/mL of mucus | FISH on Carnoy-fixed sections | directly comparable |
| CFU per biopsy | culture of a washed biopsy | divide by biopsy area, then by slab thickness |
| CFU/g | culture/qPCR of tissue or contents | assume 1 g/mL, then use the areal conversion |

For the model slab, the exact volumetric-to-areal conversion is:

```text
rho_areal [cells/cm²] = rho_vol [cells/mL] * Lz_cm
```

For the standard 100 µm slab, `Lz_cm = 0.01 cm`, so:

```text
rho_areal = rho_vol / 100
```

This factor of 100 is why an envelope quoted per cm² or per g must not be
pasted directly into `dysbiosis_threshold`.

The metric is a **domain mean over the mucus slab**. It averages local
patchiness and is therefore not meaningful as a proxy for a microcolony-scale
measurement. Patchiness should instead be compared using areal density and
patch-size distributions.

Fecal CFU/g figures do not apply: GutIBM has no lumen compartment and models
the outer loose colonic mucus layer, not fecal contents. A cells/g conversion
would additionally require an explicit mucus mass-density assumption; the
companion geometry report uses 1 g/mL only as a transparent conversion
assumption, not as a model constant.

## 2. Literature anchors

### Total mucosa-associated bacteria

- Swidsinski et al. 2005, *Journal of Clinical Microbiology* 43(7):3380–3389:
  FISH on Carnoy-fixed sections found mucosal bacteria above `10⁹/mL` in 35%
  of healthy controls, 65% of IBS, and 90–95% of IBD/self-limiting colitis.
  Mean mucosal biofilm density was approximately two logs higher in IBD than
  in controls/IBS. Carnoy fixation is load-bearing; formalin-fixed sections
  showed no biofilm.
- A 2020 *Mucosal Immunology* review, Table 1, summarizes field biofilm
  criteria: Swidsinski uses at least `0.4 × 10⁹ bacteria/mL` spanning at
  least 50 µm within 1 µm of epithelium; the Sears lab uses at least
  `2 × 10⁹/mL` spanning 150–200 µm.
- Washed healthy colon biopsies are reported to be nearly free of bacteria,
  while IBD biopsies remain dense after washing. Culture-after-washing
  therefore systematically under-recovers the healthy mucosal signal.

### Enterobacteriaceae share

- PLOS ONE 2011;6:e25042, pyrosequencing of five healthy subjects across
  colon locations, measured Proteobacteria at `8.6% ± 4.7%` of the
  mucosa-associated community.
- *Spatial Characteristics of Colonic Mucosa-Associated Gut Microbiota in
  Humans* (2021), 27 polyp-free adults across five segments, measured
  Proteobacteria at approximately 6%; Escherichia/Shigella was among genera
  above 1% relative abundance.

Enterobacteriaceae are a subset of Proteobacteria. An order-of-magnitude
healthy share of approximately 1% is therefore used as an envelope anchor:

```text
10⁹/mL total mucosal bacteria × approximately 1% ≈ 10⁷ Enterobacteriaceae/mL
```

### E. coli culture anchor

An *Inflammatory Bowel Diseases* 2013 study (30 Crohn's disease, 15
ulcerative colitis, and 14 healthy controls) measured mucosa-associated
*E. coli* using a gentamicin-protection assay on 5-mm forceps biopsies.
Healthy controls had a median 230 CFU/biopsy, compared with 1,350 for
ulcerative colitis and 2,750 for Crohn's disease. A 5-mm biopsy is
approximately 0.2 cm², giving approximately `10³ CFU/cm²`, or approximately
`10⁵ cells/mL` when represented by a 100 µm slab. Intracellular *E. coli* was
detected in 90% of Crohn's disease, 47% of ulcerative colitis, and 0% of
healthy controls.

This is approximately two logs below the FISH-derived estimate, as expected
from culture-after-washing under-recovery. The `10⁵` and `10⁷ cells/mL`
anchors are therefore treated as a healthy bracket, not competing point
estimates.

### Bloom regime

Winter & Bäumler (host nitrate, *Science* 2013, and the accompanying
*Nature Reviews Gastroenterology* commentary) describe inflammation-derived
nitrate, S-oxides, and N-oxides giving facultative anaerobic
Enterobacteriaceae a respiratory advantage over obligate anaerobes.
PLOS Pathogens 2014;10:e1003844 reports inflammation-driven blooms in which
*S. Typhimurium* and commensal *E. coli* exceed `10⁸ CFU/mL`.

That is a different ecology, not merely a denser version of the modeled
homeostatic state. Host electron acceptors, mucus erosion, and immune
infiltration are not represented here, so densities above `10⁸ cells/mL`
are extrapolation and out of scope.

### Patchiness

The reported 10–100 µm microcolony scale is consistent with biofilm criteria
spanning at least 50 µm (Swidsinski) to 150–200 µm (Sears), often within 1 µm
of epithelium. With `dx = 2–5 µm`, the model resolves a 50–200 µm patch as
approximately 10–100 cells across. The appropriate patchiness validation
target is areal density plus patch-size distribution, not the domain-mean
volumetric density.

## 3. Proposed envelope

The model's primary unit is cells/mL. The areal column is the equivalent for a
100 µm slab.

| Regime | Enterobacteriaceae (cells/mL) | Equivalent (cells/cm² at 100 µm) | Interpretation |
|---|---:|---:|---|
| **Valid** | `10⁵ – 10⁷` | `10³ – 10⁵` | Culture lower bracket through approximately 1% of a `10⁹/mL` mucosal community |
| **Caution** | `10⁷ – 10⁸` | `10⁵ – 10⁶` | Enterobacteriaceae become a major community fraction; total density enters biofilm criteria (`0.4–2 × 10⁹/mL`) |
| **Out of scope** | `>10⁸` | `>10⁶` | Bloom territory requiring host inflammatory electron acceptors and other absent ecology |

This is the proposed biological shape expressed in the model's volumetric
units; the literature anchors are primarily areal or per-mass values, so using
them without the slab conversion would shift the threshold by roughly two
orders of magnitude.

The shipped default guardrail is `1e8 cells/mL`, the top of the caution band.
The halt retains its existing strict `>` comparison and global-count/global-
volume basis. A threshold at `1e7 cells/mL` would be too trigger-happy for
the small campaign geometry described below.

## 4. Where shipped configurations sit

Initial densities from `/home/ubuntu/envelope_geometry.md` are:

| Configuration family | Geometry/count | Initial density | Position |
|---|---|---:|---|
| Stage 1/2 campaign | G3, 100 agents | `4.0e6 cells/mL` | Valid |
| Stage 3 campaign | G2, 600 agents | `1.5e6 cells/mL` | Valid |
| `single_colony` | G4, 500 agents | `5.0e6 cells/mL` | Valid |
| `eari_vadi_validation` / `fork_profile` | G3, 160 agents | `6.4e6 cells/mL` | Valid |
| `diversity_paradox` | G2, 600 agents | `1.5e6 cells/mL` | Valid |
| `cell_biology` | G1, 150 agents, 25 µm slab | `2.4e9 cells/mL` | Out of scope by this metric; explicitly opted out |
| Scaling examples | G4/G2, 1e5/1e6 agents | `1.0e9` / `2.5e9 cells/mL` | Out of scope performance runs; explicitly opted out |
| CTest scaling smoke | G5, 3,000 agents | `7.5e8 cells/mL` | Out of scope CI smoke geometry |

The out-of-scope rows are microcolony or performance configurations where a
domain-mean density is not the right biological metric. The envelope is a
property of the science configurations, not a global invariant over every
test and benchmark input.

### Small-domain guardrail consequence

In the stage-1/2 G3 geometry (`2.5e-5 mL`), the entire valid band corresponds
to only **2.5–250 agents**. The stage-1/2 configurations start at 100 agents.
A threshold at the top of the *valid* band would therefore fire during
ordinary healthy growth. The `1e8 cells/mL` threshold corresponds to 2,500
agents in G3 and 40,000 agents in G2, giving the guardrail useful dynamic
range without triggering on ordinary starting populations.

The lower end of the valid band, `1e5 cells/mL`, corresponds to only 0.025
agents in G3 and is not representable as an integer population. This is a
resolution limitation of the small-domain configuration, not a claim that
the biological lower bound is zero.

## 5. Runtime guardrail and reporting

`dysbiosis_threshold` defaults to `1e8 cells/mL`; `0` continues to mean
disabled for explicitly opted-out performance or microcolony configurations.
The campaign configurations that previously used inert `1e12` values now use
`1e8`.

The runtime progress and summary output report the current global density in
cells/mL alongside the global agent count. This uses the already-computed
global count and adds no MPI collective or HDF5 field.

The adaptive-timestep overcrowding check is intentionally separate and is not
changed by this operating-envelope update.
