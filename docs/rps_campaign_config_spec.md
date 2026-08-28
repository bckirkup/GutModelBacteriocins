# Rock-Paper-Scissors Campaign Config Spec

# Literature-grounded parameterization for E. coli colicin dynamics

## Parameter Summary (from literature review)

|Parameter|Best estimate|Range|Evidence|Citation|
|-|-|-|-|-|
|BtuB-null growth cost|2% mu\_max reduction|0-10%|Inferred; condition-dependent|Feldgarden \& Riley|
|FepA-null growth cost|5%|0-20%|Inferred; iron-dependent|Feldgarden \& Riley|
|Double-null (BtuB+FepA)|8%|0-30%|Approx multiplicative|--|
|Immunity protein cost|0.2%|0-2%|Inferred; below assay resolution|ColIb study|
|Fecal total corrinoids|0.97 uM (0.97 mmol/m3)|0.2-2.4 uM|Direct measurement|Allen \& Stabler|
|True cobalamin|14 nM (0.014 mmol/m3)|3-100 nM|Direct; 1.4% of total|Allen \& Stabler|
|Mucus corrinoid conc|0.3 uM total|0.03-1 uM|Inferred; partition \~0.3|No direct measurement|
|BtuB-cobalamin Kd|0.3 nM|0.05-5 nM|Literature-based|Cadieux, structural studies|
|BtuB-colicin E Kd|1 nM|0.1-20 nM|Literature-based|Structural/thermodynamic|
|Analogue-BtuB Kd|10 nM|0.3 nM - 1 uM|Sparse data|--|
|Immigration rate|0.005/day per patch|1e-4 to 0.1/day|Inferred from exposure \~5 strains/host/day|Caugant, Gordon|
|Concurrent strains|2 dominant|1-4|Historical; undersampled|Caugant 1981|
|Mucus E. coli density|1e7 cells/g mucus|1e4-1e9|Inferred; no direct human measurement|--|
|Areal density|1e5 cells/cm2|1e2-1e7|Derived from volume density|--|
|Colicin D\_free (water)|8.5e-11 m2/s|7-10e-11|Analogue (BSA/ovalbumin)|--|
|Colicin D\_mucus|7e-11 m2/s|2e-11 - 1.2e-10|Analogue (cervical mucus)|Cone 2009|
|Mucus retardation factor|1.2 (D\_mucus/D\_water \~ 0.68-1.0)|1-5|Cervical mucus analogue|--|
|SOS induction rate|1% per generation|0.5-3%|Reporter assays|ColIb reporter|
|Lysis fraction of induced|20%|5-100%|System-dependent|ColIb/prophage study|
|Induction-to-lysis delay|90 min|30-240 min|Inferred|--|
|Burst size (molecules)|1e5|1e3-1e6|Major data gap|No direct measurement|
|Lethal hit threshold|1 translocated molecule|1-10|Classical one-hit kinetics|--|
|Binding to commitment|5 min|1-20 min|Adsorption experiments|de Zwaig \& Luria|
|ColE1 self-transfer|0 (non-conjugative)|exactly 0|Established|--|
|Conjugative Col plasmid|1e-10 mL/cell/h|1e-12 to 1e-9|Sparse, highly variable|--|



## Three-Strain Config for Devin

```json
{
  "\_comment": \[
    "Rock-Paper-Scissors: Producer (P) kills Sensitive (S),",
    "Resistant (R) outgrows Producer (no lysis cost),",
    "Sensitive (S) outgrows Resistant (no resistance cost).",
    "Parameters from literature review (task:37a7d949)."
  ],
  "initial\_strains": \[
    {
      "type": 1,
      "count": 20,
      "mu\_max": 0.00055,
      "label": "Producer",
      "plasmids": \[
        "ColE1",
        "ColB"
      ],
      "conjugative": false,
      "\_note": "SOS lysis + colicin release. Immunity to own colicins. mu\_max same as Sensitive (cost is lysis, not growth)."
    },
    {
      "type": 2,
      "count": 20,
      "mu\_max": 0.00055,
      "label": "Sensitive",
      "plasmids": \[],
      "conjugative": false,
      "\_note": "No immunity. Killed by colicin. Fastest grower (no costs)."
    },
    {
      "type": 3,
      "count": 20,
      "mu\_max": 0.000539,
      "label": "Resistant",
      "plasmids": \[],
      "conjugative": false,
      "receptor\_mutations": {
        "BtuB": "null",
        "FepA": "null"
      },
      "\_note": "BtuB+FepA null = immune to ColE + ColB. 2% growth cost = mu\_max \* 0.98. Literature: 0-10% for BtuB alone."
    }
  ],
  "receptor\_config": {
    "kd\_colicinE\_btuB": 1e-09,
    "kd\_b12\_btuB": 0.0001,
    "kd\_colicinB\_fepA": 2e-09,
    "kd\_enterobactin": 1e-06,
    "kill\_rate\_colicin": 0.001,
    "\_note": "kd\_b12\_btuB = 1e-4 gives 11x competitive factor at 1e-3 mol/m3 corrinoid. From corrinoid sweep: graded kills, not binary."
  },
  "corrinoid\_field": {
    "concentration\_mol\_m3": 0.0003,
    "\_note": "0.3 uM total corrinoids in mucus. Inferred from fecal 0.97 uM \* partition 0.3. True Cbl is \~14 nM but analogues dominate."
  },
  "immigration": {
    "rate\_per\_day": 0.005,
    "types\_introduced": \[
      1,
      2,
      3
    ],
    "count\_per\_event": 1,
    "\_note": "\~5 strains/host/day exposure, but establishment prob \~0.001. Inject one cell of random type every \~200 sim-days per patch."
  },
  "colicin\_parameters": {
    "sos\_lysis\_prob": 0.01,
    "burst\_size\_molecules": 100000.0,
    "burst\_release\_tau\_s": 300,
    "D\_free\_colicin\_m2s": 7e-11,
    "retardation\_mucus": 1.2,
    "toxin\_cutoff\_m": 0.0002,
    "\_note": "D\_free from cervical mucus analogue. Retardation 1.2 (not 10). 10x was too high per Cone 2009: Dmucus/Dwater \~ 0.68-1.0 for 40-70 kDa proteins."
  },
  "delivery": {
    "mode": "flux",
    "epithelial\_flux": 3.2268e-09,
    "\_note": "flux\_0.3x from delivery arms. Runs 24h inside population bound."
  },
  "vbf": {
    "carbon\_sink\_vmax": 5.5e-05,
    "\_note": "1x default. Binary switch (on/off), not a dial."
  }
}
```

## Campaign Design

### Core sweep: resistance cost x corrinoid competition

|Factor|Values|Rationale|
|-|-|-|
|Resistant mu\_max cost|0%, 2%, 4%, 8%|Literature range 0-10%; need to find cycling boundary|
|kd\_b12\_btuB|1e-5, 1e-4, 1e-3|Competition factor 101/11/2; controls colicin efficacy|
|Colicin retardation|1.2, 3.0, 10.0|Literature supports 1.2; 10 was previous default|

4 cost x 3 kd x 3 retardation = 36 cells
x 5 seeds x 7-day horizon = 180 runs

### Immigration panel

At best RPS operating point from core sweep:
Immigration rate: \[0, 0.001, 0.005, 0.02, 0.1] /day
5 rates x 5 seeds = 25 runs

### Key outputs

1. Does cycling occur? (all 3 types persist for >48h)
2. At what resistance cost does P displace R?
3. At what kd does colicin become irrelevant?
4. What immigration rate sustains 3-type coexistence?
5. Kills/lysis ratio across the parameter space

### Critical parameter change from current defaults

|Parameter|Current|Literature|Change|
|-|-|-|-|
|Colicin retardation|10.0|1.2 (Dmucus/Dwater \~ 0.8)|**Reduce 8x**|
|kd\_b12\_btuB|1e-6|1e-4 (11x competition at 0.3 uM corrinoid)|**Increase 100x**|
|Corrinoid field|1e-3 mol/m3|3e-4 mol/m3 (0.3 uM mucus estimate)|**Reduce 3x**|
|Resistant strain mu\_max|not present|5.39e-4 (2% cost)|**New strain type**|

The retardation change is the biggest: from 10 to 1.2.
This increases the colicin kill radius from \~22 um to \~60 um per step,
and the steady-state killing zone by \~7x in area.
Combined with kd\_b12\_btuB = 1e-4 (11x competition instead of 1001x),
this should make colicin a functional weapon rather than decorative.

