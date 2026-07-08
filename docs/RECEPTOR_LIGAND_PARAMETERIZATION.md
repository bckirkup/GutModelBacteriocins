# Receptor Ligand Parameterization

This document grounds the competitive-binding parameters used by
`FixReceptor` in gut biochemistry. It summarizes the literature review that
motivated the Spec 6 nutrient-cycle rework and records the values adopted in
the model, along with the key unknowns left for future work.

## Competitive-binding model

Colicin kill probability at each receptor follows competitive Michaelis–Menten
binding between the toxin and the receptor's native ligand:

```
apparent_Kd = kd_tox × (1 + [ligand] / kd_ligand)
occupancy   = receptor_expr × [tox] / (apparent_Kd + [tox])
P(kill)     = 1 − exp(−kill_rate × occupancy × immunity_eff × dt)
```

Implemented in `FixReceptor::toxin_occupancy` (`src/fixes/fix_receptor.cpp`).

## 1. BtuB — corrinoids vs the colicin E family

The gut contains ~1309 ng/g total corrinoids (Allen & Stabler 2008), but true
cobalamin is only ~1.4% of the pool. The dominant species are non-cobalamin
analogs ([2-methyladenine]cobamide 60.6%, [p-cresol]cobamide 16.3%,
pseudo-B12 12.5%). Converting to molarity (mean cobamide MW ~1350 Da):

- Total corrinoids: ≈ **1 µM** (9.7e-7 mol/m³)
- True cobalamin only: ≈ 14 nM
- Dominant analog: ≈ 590 nM

Because BtuB is both the corrinoid importer and the colicin-E receptor, ambient
corrinoid competitively protects cells from colicin E. Functional competition
is confirmed: 5 µM B12 preincubation blocks ColE1 killing (Masi et al. 2007).

| Parameter | Model value | Basis |
|-----------|-------------|-------|
| `[B12]` field | **1e-6 mol/m³ (1 µM)** | Total bioavailable corrinoid pool, not true Cbl only |
| `kd_b12_btuB` / `kd_corrinoid_btuB` | 1e-9 (1 nM) | Sub-nM transport affinity supported; **key unknown** for the dominant analog |
| `kd_colicinE_btuB` | 5e-10 (0.5 nM) | No better data; colicins operate at pM–nM potency |

**Consequence.** At `[corrinoid]` = 1 µM and `kd_ligand` = 1 nM the competitive
factor is `1 + 1e-6/1e-9 ≈ 1001`, so `apparent_Kd ≈ 0.5 µM`. Colicin E is ~1000×
less potent than under the old `[B12]` = 1 nM field (competitive factor ≈ 2).
This is the single most consequential downstream effect of the nutrient rework.

**Key unknown / future sweep.** BtuB affinity for the dominant analog
([2-MeAde]cobamide) is unmeasured; if it binds 10×/100× weaker, the competitive
factor falls to ~101×/~11×. A sweep over
`kd_corrinoid_btuB ∈ {1e-9, 1e-8, 1e-7, 1e-6}` at `[corrinoid]` = 1 µM is
recommended future work (out of scope for the Spec 6 PR).

## 2. FepA — siderophores vs colicin B/D

Live-cell FepA–ferric-enterobactin adsorption Kd is ~0.2 nM (Klebba 2003);
colicins B/D share FeEnt binding determinants at extracellular loop 5 (true
competitive inhibition). Under iron limitation the Fur regulon is induced 35–56×
(cirA 56×, fepA 35×, fiu 39×), though those figures come from antimicrobial-
peptide stress rather than pure iron chelation.

| Parameter | Model value | Basis |
|-----------|-------------|-------|
| `kd_enterobactin` | **1e-9 (1 nM)** | Live-cell Kd ~0.2 nM; 1 nM is conservative (was 10 nM) |
| `kd_colicinB_fepA` | 2e-9 (2 nM) | No data; same order as colicin E at BtuB |
| `fur.upregulation_max` | **10.0** | Conservative vs measured 35–56×; was 4.0 |

## 3. CirA — linearized enterobactin vs colicin Ia / microcin V

No quantitative Kd data for CirA were found; current estimates are retained.
`cirA_linearized_fraction` (default 0.3) remains a tunable approximation of the
esterase-linearized enterobactin fraction. Note the natural "Trojan horse"
salmochelin–microcin H47 conjugate (Cherrak 2024) that the per-receptor toxin
field system could capture if IroN-targeted microcins are added to the plasmid
library.

## 4. Key unknowns and future experimental targets

| Unknown | Impact | How to resolve |
|---------|--------|----------------|
| BtuB Kd for [2-MeAde]cobamide | Determines colicin E effectiveness (1000× range) | ITC/SPR: purified BtuB + gut corrinoid analogs |
| FepA Kd for colicin B | Colicin B effectiveness | SPR: FepA + ColB |
| CirA Kd for colicin Ia | Colicin Ia effectiveness | SPR: CirA + ColIa |
| Gut free-iron concentration | Competitive binding baseline | ICP-MS on mucosal samples |
| Corrinoid spatial gradient | Position-dependent protection | LC-MS on fractionated gut samples |

## References

- Allen & Stabler 2008. Am J Clin Nutr 87:1324–1335.
- Masi et al. 2007. (B12 competition with ColE1 at BtuB.)
- Mok et al. 2025. Laboratory evolution of E. coli with pseudocobalamin.
- Klebba 2003. Front Biosci 8:s1422–36.
- Thulasiraman et al. 1998. J Bacteriol 180:6689–96.
- Cherrak et al. 2024. Nat Microbiol. doi:10.1038/s41564-024-01684-1.
- Xiao et al. 2017. Annu Rev Nutr 37:103–130.
- Degnan et al. 2014. Cell Metab 20:769–778.
