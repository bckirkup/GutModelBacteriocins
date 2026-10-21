# GutIBM: Reconciliation of Molecular vs Culture Mucosal E. coli Densities

## The contradiction

Two independent datasets give mucosal Enterobacteriaceae/E. coli densities
differing by ~100–1000×:

| Source | Method | Value | Converts to |
|---|---|---|---|
| Ahmed et al. 2007 | qPCR (16S copies/mg tissue) | 10^5.5 Enterobacteriaceae copies/mg | ~4.5×10⁷ cells/mL tissue (÷7 rRNA operons) |
| Elliott et al. 2013 | Selective culture (CFU/biopsy) | 230 CFU/20mg biopsy = 1.15×10⁴ CFU/g | ~4.6×10⁴–2.3×10⁵ CFU/mL mucus |

Gap: ~100–1000×.

## Why both are correct

This gap is well-documented across gut microbiology and reflects real
differences in what each method measures:

### 1. Taxonomic resolution
- **qPCR**: Targets the Enterobacteriaceae *family* (16S primers). Includes
  Klebsiella, Enterobacter, Citrobacter, Proteus, Escherichia, plus
  unculturable lineages. E. coli is a subset.
- **Culture**: MacConkey agar selects for viable E. coli specifically (or at
  most Enterobacteriaceae that grow on selective media).

### 2. Viability
- **qPCR**: Detects DNA from viable cells, dead cells, VBNC (viable but
  non-culturable) cells, free extracellular DNA, and phage-associated DNA.
  No viability discrimination.
- **Culture**: Counts only cells that can form colonies on selective media.
  Cells in stationary phase, biofilm state, or nutrient-stressed VBNC state
  (common in mucus: 10–100× lower culturability) are missed.

### 3. Compartment
- **qPCR**: Per mg of biopsy *tissue* — includes epithelium, lamina propria,
  submucosal tissue, plus whatever mucus and bacteria are retained after
  processing. The 16S copies are per mg of this mixed material.
- **Culture**: Mucosa-*associated* bacteria from washed biopsies — surface
  and loosely adherent organisms. More specifically mucosal.

### 4. Copy number correction
- Enterobacteriaceae carry ~7 rRNA gene copies per genome (range 5–8).
  Raw 16S copy counts must be divided by ~7 to estimate cell counts.
  This accounts for ~1 log but not the full 2–3 log gap.

## What each dataset tells the model

| Dataset | What it constrains | What it does NOT tell us |
|---|---|---|
| Ahmed qPCR | **Relative fraction** of Enterobacteriaceae along the colon (~1–2%, no significant proximal-distal gradient, P=0.09) | Absolute viable E. coli density in mucus |
| Elliott culture | **Absolute viable E. coli density** in healthy vs inflamed mucosa (healthy ~230 CFU/biopsy, CD ~2750 CFU/biopsy) | Whether Enterobacteriaceae fraction varies by segment |

The two datasets are **complementary, not contradictory**. Ahmed constrains
the model's spatial uniformity assumption (no preferential proximal
colonization). Elliott constrains the model's absolute density scale and
the healthy-to-inflamed range.

## Resolution for the model

The GutIBM simulates **viable, metabolically active E. coli** in outer mucus.
This corresponds to the culture-based measurement.

### Operating range (culture-anchored)

| Regime | Viable E. coli (CFU/mL mucus) | Source |
|---|---|---|
| Healthy | 10⁴–10⁵ | Elliott healthy biopsies |
| Transition/alert | 10⁵–10⁶ | Elliott UC range |
| Dysbiosis guard | 10⁶ | Elliott inflamed CD |
| Model-invalid | ≥10⁷ | Beyond any healthy/early-inflammatory culture value |

### Enterobacteriaceae fraction (molecular-anchored)

| Constraint | Value | Source |
|---|---|---|
| Fraction of mucosal community | ~1–2% | Ahmed qPCR |
| Proximal-distal gradient | Not significant (P=0.09) | Ahmed qPCR |
| Individual variation | 0.1–5% | Ahmed qPCR (SD spans order of magnitude) |

### Stool target (culture-based)

| Quantity | Value | Basis |
|---|---|---|
| E. coli in healthy stool | 10⁶–10⁸ CFU/g | Selective plating, standard reference |
| Total fecal bacteria | ~10¹¹ cells/g | Molecular (FISH/flow cytometry) |
| E. coli fraction of stool | ~0.001–0.1% | Ratio of above |

## The molecular number in context

The ~4.5×10⁷ Enterobacteriaceae cells/mL tissue from qPCR should be
understood as:

- Includes all family members, not just E. coli
- Includes non-viable cells and free DNA
- Is per mg of tissue (mixed epithelial + stromal + mucosal material)
- Is ~100× higher than the viable E. coli count in the same environment

This is consistent with the general finding across gut microbiology that
molecular methods detect 10–100× more bacteria than culture
(Eckburg et al. 2005, Suau et al. 1999).

For the model: the qPCR number defines the *total Enterobacteriaceae signal*
that the VBF implicitly includes. The culture number defines the *viable
E. coli* that the model explicitly simulates. The model's VBF density
(currently 10¹¹ cells/m³ = 10⁸ cells/mL) represents the molecular-scale
total community, not the viable-culture-scale community.

## Consequences for Spec 13

1. **Layer 3 regional table**: Use Ahmed fractions (no gradient) but do NOT
   convert Ahmed absolute 16S copies to model cell densities. The absolute
   scale comes from Elliott.

2. **Dysbiosis guard**: Anchored to culture at 10⁶ CFU/mL viable E. coli.
   This is within ~5× of Elliott inflamed CD mucus values.

3. **Stool prediction**: The wall-to-stool equation predicts viable CFU/g,
   validated against culture-based stool counts (10⁶–10⁸ CFU/g), not against
   molecular totals.

4. **VBF density**: The existing 10¹¹ cells/m³ (10⁸ cells/mL) is consistent
   with FISH-based total mucosal bacteria (~2×10⁸ cells/mL in healthy
   controls, Swidsinski 2005). This is the right scale for the VBF.

## Layer 3 fraction denominator — resolved 2026-10-21

The Layer-3 Enterobacteriaceae fraction readout
(`enterobacteriaceae_fraction = segment_mean_cfu_ml / flora_density_cfu_ml`)
shipped with `flora_density_cfu_ml = 1e7`, which reported 0.14–0.41
mucosal fractions — against a ~1–3% empirical band that suggested a real
physics defect. It was a denominator error: the proxy had been set an
order of magnitude below the VBF/FISH scale this doc already derives.

Literature check (Consensus deep-research thread, record committed at
`analysis/spec13_gradient/literature_mucosal_density_consensus.md`):

- **Absolute cells/mL of native human colonic mucus is not directly
  measured.** Closest measured quantities: Earley 2020 copies/mg mucus
  (RT-qPCR), Poxton 1997 CFU/biopsy (aerobes 2.4e3–1.3e6, anaerobes
  1.4e5–3e7 per ~5.6 mg biopsy), Zoetendal 2002 1e5–1e6 cells/biopsy
  (flow cytometry), Lyra 2012 log10 genomes/g wet tissue (qPCR),
  Swidsinski FISH ~0.02×10^10/mL ≈ 2e8 cells/mL healthy controls.
- **Mucosal Enterobacteriaceae fraction**: Hong 2011 (16S
  pyrosequencing, Table 2) — 0.68–3.81% across 5 healthy individuals;
  consistent with Ahmed's 1–2%. Proteobacteria fraction: 8.6% (Hong),
  11% (Shen), 19.3% distal mucosa vs 9.3% feces (Ringel).

Resolution: `flora_density_cfu_ml` is a **FISH-scale total-community
cells/mL** denominator, shipped at 1e8 (declared interval 1e8–2e8,
Swidsinski + VBF scale). At that denominator the audited spec13 runs'
0.14–0.41 fractions become 1.4–4.1% (0.7–2% at 2e8) — inside or at the
top of the Hong/Ahmed empirical band. The fraction *level* was a
measurement-scale defect in the readout denominator, not a physics
error; `patch_supply_scale` needed no correction for this question.

## Generalization

This is not specific to density. The same substitution — a precise
sequence-derived number standing in for an imprecise functional one — recurred
in the Spec 14 lysogen prior, where integrated-prophage detection in MAGs was
adopted as the prevalence of *inducible* prophage. See
`docs/SPEC14_PRIOR_REVIEW.md` §4 for the class and the rule.

## References

Ahmed et al. 2007. Mucosa-associated bacterial diversity in relation to human
terminal ileum and colonic biopsy samples. Appl Environ Microbiol.
doi:10.1128/aem.01143-07

Elliott et al. 2013. Quantification and characterization of mucosa-associated
and intracellular E. coli in IBD. Inflamm Bowel Dis.
doi:10.1097/mib.0b013e3182a38a92

Swidsinski et al. 2005. Spatial organization of bacterial flora in normal and
inflamed intestine. J Clin Microbiol.
doi:10.1128/jcm.43.7.3380-3389.2005
