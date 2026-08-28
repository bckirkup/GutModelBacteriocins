# Variant Surveillance System Spec
# CTB Paper 3: Cruise Ships as Phylogenomic Observatories

## Overview

Extend CTB with heritable strain identity, within-voyage evolution, and
sequencing observation models to study how intensive shipboard monitoring
can detect novel pathogen variants and characterize transmission parameters
faster than land-based surveillance.

Three capabilities:
1. **Strain tracking**: heritable lineage labels with mutation at transmission
2. **Variant observation**: amplicon sequencing, Nanopore, wastewater deconvolution
3. **Surveillance economics**: cost-benefit of onboard vs ashore investment

## 1. Heritable Strain System

### 1.1 Strain identity on infection

Every infection event carries a StrainState with:
- strain_id: unique lineage identifier
- parent_strain_id: who transmitted to me
- generation: transmission generations from founder
- transmissibility_multiplier: dose_adj modifier (1.0 = reference)
- shedding_multiplier: peak shedding amplitude modifier
- incubation_modifier: days added/subtracted from reference
- immune_escape: 0.0 = full cross-immunity, 1.0 = no cross-immunity
- n_mutations: accumulated from founder
- genotype: e.g. GII.4 vs GII.17 for norovirus

### 1.2 Mutation at transmission

When agent A infects agent B, the child strain inherits the parent with:
- Bernoulli(mutation_rate) chance of a point mutation
- If mutated: Bernoulli(phenotype_mutation_fraction) chance it affects phenotype
- Phenotype targets: transmissibility, shedding, incubation, immune_escape
- Recombination deferred to v2 (requires co-infection mechanics)

### 1.3 Per-pathogen strain parameters

| Pathogen | Mutation rate/transmission | Phenotype fraction | Genotypes |
|---|---|---|---|
| Norovirus GII | 0.02 | 0.05 | GII.4, GII.17, GII.2 |
| SARS-CoV-2 | 0.03 | 0.10 | by lineage |
| Influenza A | 0.04 | 0.15 | H1N1, H3N2 |
| Measles | 0.005 | 0.01 | genotypes A-H |
| Rigelian Fever | 0.05 | 0.20 | Alpha, Beta, Gamma |
| Psi-2000 Polywater | 0.10 | 0.50 | Phase I, Phase II |
| Barclay Protomorphosis | 0.05 | 0.30 | Genesis-prime, Genesis-delta |
| TNG Shipboard Influenza | 0.04 | 0.15 | -- |

Trek pathogens are parameterized MORE mutationally dynamic than real ones.
Psi-2000 especially (rapid behavioral modification = extreme adaptability).

### 1.4 Cross-immunity matrix

Per pathogen: NxN matrix where entry [i,j] = protection from genotype i against j.
Norovirus: same-genotype ~0.85, cross-genotype ~0.15-0.20.
Immune escape mutations reduce effective protection: base * (1 - immune_escape).

## 2. Sequencing Observation Models

### 2.1 Clinical specimen sequencing

| Pathogen | Amplicon target | Nanopore accuracy | Ct threshold | Cost |
|---|---|---|---|---|
| Norovirus | Capsid (VP1) | 92% | <30 | $150 |
| SARS-CoV-2 | Spike + WGS | 97% | <28 | $100 |
| Influenza | HA + NA | 95% | <30 | $120 |
| Measles | N gene | 98% | <25 | $100 |
| Rigelian Fever | 16S + virulence | 90% | <32 | $200 |
| Psi-2000 | BML (behavioral mod locus) | 85% | <20 | $300 |
| Barclay Protomorphosis | Genesis factor | 88% | <25 | $250 |

### 2.2 Wastewater strain deconvolution

Freyja-style: wastewater contains a mixture of strains from all shedders.
True proportions weighted by shedding rate.
Observation: Dirichlet-multinomial draws scaled by sequencing depth.

### 2.3 Surface sampling strain recovery

Surface contamination carries strain_id of depositing agent.
Recovery probability depends on surface type, time since deposition.

## 3. Federation Port Network

| Port | Region | Pop | Surveillance | Endemic pathogens | WBE | Analog |
|---|---|---|---|---|---|---|
| Starbase 1 (Earth) | Sol | 5000 | Full | Influenza analog | Yes | Miami |
| Starbase 11 | Alpha Core | 2000 | Standard | Rigelian Fever | Yes | Barcelona |
| Farpoint Station | Alpha Frontier | 800 | Limited | Protomorphosis | No | Roatan |
| Deep Space 9 | A-G Border | 3000 | Mixed | Rigelian, Psi-2000 | Yes | Singapore |
| Starbase 74 | Alpha | 1500 | Standard | -- | Yes | Copenhagen |
| Risa | Alpha | 50000 | Minimal | Influenza, Psi-2000 | No | Cozumel |
| Memory Alpha | Alpha | 200 | None | -- | No | Skagway |
| K-7 Station | Klingon Border | 500 | Limited | Rigelian Fever | No | Juneau |

Risa is the Federation's Cozumel: high mixing, poor surveillance, tourism.
DS9 is the epidemiological chokepoint (wormhole = novel pathogen introduction).
Memory Alpha is Skagway (tiny, no infrastructure).

### Starfleet itineraries

Constitution (TOS): Starbase 1 -> Starbase 11 -> K-7 -> Memory Alpha -> Starbase 1
Galaxy (TNG): Starbase 1 -> Farpoint -> DS9 -> Risa -> Starbase 74 -> Starbase 1

## 4. Surveillance Economics

| Scenario | Onboard | Ashore | Annual cost | Variants detected |
|---|---|---|---|---|
| Baseline | Syndromic only | VSP reporting | ~$0 | 0 |
| Minimal | Clinical qPCR | Reference lab | ~$50K/ship/yr | 2-5 |
| Moderate | Nanopore + WW qPCR | Regional lab | ~$150K/ship/yr | 5-15 |
| Full | Nanopore + WW amplicon + surface | Bioinfo hub + port WBE | ~$500K/ship/yr | 15-30 |
| Fleet network | Full x 10 ships | Centralized + port network | ~$3M/yr | 50-100 |

## 5. Implementation Plan

Phase 1 (PRs 1-3): Strain tracking — StrainState, mutation, cross-immunity
Phase 2 (PRs 4-6): Observation models — amplicon, Nanopore, WW deconvolution
Phase 3 (PRs 7-9): Port networks + economics — Federation ports, cost model, campaigns

## 6. Campaigns

| Campaign | Runs | Purpose |
|---|---|---|
| Variant emergence (real) | 2,000 | Detect novel norovirus/SARS-CoV-2 variants |
| Federation patrol (Trek) | 1,000 | DS9 wormhole variant introduction |
| Fleet network value | 3,000 | 10-ship fleet vs single ship detection power |
| Investment optimization | 2,000 | Cost-benefit frontier onboard vs ashore |
| Risa outbreak | 500 | Tourism planet, poor surveillance |

## 7. Paper 3 Outline

Title: Cruise Ships as Phylogenomic Observatories

1. Introduction: the case for ship-based variant surveillance
2. Methods: strain tracking, sequencing observations, cost framework
3. Results:
   a. Within-voyage evolution: norovirus genotype dynamics in 7-day outbreaks
   b. Detection speed: ship vs land across pathogens
   c. Transmission parameter estimation from closed cohorts
   d. Fleet network detection power
   e. Investment frontier
4. Generalization: Federation port networks
5. Discussion: drift, bottlenecks, ships as natural experiments
6. Appendix: Starfleet variant surveillance
