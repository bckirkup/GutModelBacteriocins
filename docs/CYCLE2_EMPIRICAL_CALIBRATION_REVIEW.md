# Cycle 2 Empirical Calibration Review (Consensus Literature Synthesis)

**Status:** Draft for review and revision  
**Target Repository:** `bckirkup/GutModelBacteriocins`  
**Companion Audit:** [SPEC13_14_MODEL_AUDIT.md](SPEC13_14_MODEL_AUDIT.md)  
**Companion Data Tables:**
- [consensus_empirical_benchmarks.csv](consensus_empirical_benchmarks.csv) (41 harmonized empirical benchmarks across 4 pillars)
- [parameter_provenance_audit.csv](parameter_provenance_audit.csv) (23 parameter-level provenance entries)
- [open_gaps_and_validation_matrix.csv](open_gaps_and_validation_matrix.csv) (12 gaps, `GAP-01..GAP-12`)

---

## 1. Executive Summary

This review synthesizes empirical literature benchmarks queried via the Consensus API (210 unique DOIs indexed; 24 full texts, 73 full-text-chunk records, 101 abstract-only records, and 12 publisher-blocked full-text checks) to resolve four Cycle 1/Cycle 2 calibration gaps in GutIBM:

1. **Free vs. Total Colonic Vitamin B12 (`GAP-01`, `b12.initial_conc`):** Shipped `1.0e-3 mol/m^3` (`1 uM`) matches **total** human fecal corrinoids (`9.66e-4 to 2.00e-3 mol/m^3`), which are predominantly non-cobalamin analogues (`>98%`) and intracellular microbial stores. Measured **true cobalamin** in human feces is `1.35e-5 mol/m^3` (`1.4%` of total), and baseline mouse luminal/fecal active cobalamin spans `4.43e-8 to 7.38e-6 mol/m^3`—74x to 22,600x lower than the shipped default. Treating `1.0e-3 mol/m^3` as free extracellular CNCbl competitively saturates BtuB (`kd_b12_btuB = 1.0e-6 mol/m^3`) and blocks ColE2 exclusion in Spec 14 V1.
2. **Mucosal Oxygen & Anaerobic Maintenance (`GAP-08`, `chem_env.oxygen.*`):** Shipped `epithelial_conc = 55.0e-6 mol/m^3` (`~0.042 mmHg`) sits 71x–636x below mammalian colonic mucosal measurements (`3–11 mmHg` / `5–35 uM` = `3.93e-3 to 3.50e-2 mol/m^3`). Conversely, `anaerobic_maintenance_factor = 15.0` (which scales `carbon_maintenance_rate` in `src/fixes/fix_metabolism.cpp:662-667`) is 2.6x–4.4x too high: *E. coli* K-12 growth-independent ATP maintenance is slightly *lower* anaerobically (`0.054` vs. `0.075 mol ATP/C-mol/h`, ratio `0.72`), and with submaximal aerobic `P/O = 1.0–1.5` (`~14–20` vs. `~2.5–3.0` anaerobic mixed-acid ATP/Glc), the empirical carbon-maintenance factor is **`~3.4–5.8x`** (chemostat `q_Glc` ratio `4.61x`).
3. **Host-Specific GI Transit & Colicinogeny Scaling (Pillar 3, `layer3.transit_h`):** Shipped `layer3.transit_h = 12/8/10 h` (`30.0 h` total) matches the upper boundary of modern human ingestible-capsule colonic transit medians (`13.6–28.9 h`) and growing-pig colonic passage (`16.5–48.8 h`), but is 6x–10x slower than murine colonic transit (`~1–5 h`; WGTT `2–8 h`). Per Barnes, Sidhu & Gordon (2007), fast gut turnover suppresses colicin producers unless mucus/crypt retention preserves local toxin concentration.
4. **ColE2 Exclusion & Prophage Induction Ceilings (Pillar 4, `GAP-10`):** Pure ColE2 exclusion fails without low resource overlap (Bakkeren et al., 2025) and nutrient-gated CsrA lysis execution (`1%` in M9-Glc vs. `15–44%` in LB; Hol et al., 2014). Proposed Spec 14 priors (`0.8 * 1.75e-7 /s = 0.0120–0.0121 /bact/day`) exceed the human gut community phage-induction ceiling (`0.001–0.01 events/bacterium/day` = `1.16e-8 to 1.16e-7 /s`, verified via `10.1016/j.cub.2025.01.006` corroborating blocked `10.1016/j.cub.2025.03.073`) before any SOS damage induction is added.

---

## 2. GAP-01: Colonic Vitamin B12 and Corrinoid Speciation

**Unit harmonization:** Cyanocobalamin `MW = 1355.37 g/mol` with luminal density `rho = 1.0 g/mL` gives:
- `1 ng/g ~ 1 ng/mL = 7.378e-7 mol/m^3` (`0.7378 nM`)
- `1 ppm (1 ug/g) = 7.378e-4 mol/m^3` (`0.7378 uM`)

| ID | Paper (Year, DOI) | Depth | System / Assay | Raw Measurement | GutIBM Value | Discrepancy vs. Baseline |
| --- | --- | --- | --- | --- | --- | --- |
| B12-01 | Robert H. Allen & Sally P. Stabler (2008, `10.1093/ajcn/87.5.1324`) | Blocked-In-FullText | Human adult feces (n=18, B12 intake <=25 ug/d) | 1309 ng cobalamin eq/g wet wt | 9.66e-04 (mol/m^3) | Matches shipped 1.0e-3 mol/m^3 (0.97x) only as total corrinoid pool |
| B12-02 | Robert H. Allen & Sally P. Stabler (2008, `10.1093/ajcn/87.5.1324`) | Blocked-In-FullText | Human adult feces (n=18, B12 intake <=25 ug/d) | 18.33 ng/g wet wt | 1.35e-05 (mol/m^3) | 74x lower than shipped 1.0e-3 mol/m^3 (0.0135x); cobinamide 1.8% |
| B12-03 | Caleb James Kelly et al. (2019, `10.1080/19490976.2019.1597667`) | Blocked-In-FullText | C57BL/6 mice (cecal vs fecal) | 0.08 cecal; 0.06 fecal ng/mg cecal; ng/mL fecal | 5.90e-05 cecal; 4.43e-08 fecal (mol/m^3) | Fecal baseline ~22,600x below shipped 1.0e-3; high oral supp (3.94 ug/mL water) gives 7.81 ng/mL fecal (5.76e-6 mol/m^3) |
| B12-04 | Andrew J. Forgie et al. (2023, `10.1186/s40168-023-01461-w`) | Full-Text | C3H/HeOuJ mice (naive & C. rodentium) | 0.01 baseline; ~15 supp ppm (ug/g) | 7.38e-06 base; 0.0111 supp (mol/m^3) | Baseline ~135x below shipped 1.0e-3 mol/m^3; ~1000x supra-physiological supp needed to reach ~1.11e-2 mol/m^3 |
| B12-05 | Palni Kundra et al. (2022, `10.3389/fnut.2022.1070155`) | Full-Text | Healthy human donors (n=8) & MacFarlane batch ferm. | 2716 median (684-4222); 5 medium ng/g feces; ng/mL medium | 2.00e-03 med; 3.69e-06 medium (mol/m^3) | Total fecal pool ~2.00e-3 mol/m^3, but standard colonic medium starts at 3.69e-6 mol/m^3 (271x below shipped) |
| B12-06 | Zachary F. Hallberg et al. (2022, `10.1021/acs.biochem.2c00367`) | Full-Text | Human fecal (n=18 avg) & bovine rumen/GI | >=7 cobamides; ~0.3 Simpson count; index | >=7 cobamide species; [2-MeAde]Cba primary (speciation note) | BtuB binds diverse cobamides/corrinoids; total corrinoid pool cannot be treated as pure free CNCbl |

**Calibration conclusion:** Kundra et al. (2022) further report *in vitro* MacFarlane fermentations without exogenous B12 (`NB12`) producing `152 +/- 4 to 564 +/- 25 ng/mL` total B12 (`Control`: `122 +/- 4 to 404 +/- 57 ng/mL`), partitioned heavily into intracellular biomass (`83 +/- 9 to 201 +/- 15 ng/mL` in D3/D5–D8 controls). Baseline free extracellular cobalamin should be parameterized at **`1.0e-6 to 1.5e-5 mol/m^3`** (`1–15 nM`; or `3.69e-6 mol/m^3` for MacFarlane `5 ng/mL`), reserving `1.0e-3 mol/m^3` strictly for total-corrinoid or bound-pool sensitivity sweeps.

---

## 3. GAP-08: Mucosal Oxygen Boundaries & Anaerobic Maintenance Scaling

### 3.1 Mucosal and Luminal `pO2` (`chem_env.oxygen.epithelial_conc`)

**Unit harmonization:** Consistent with `src/io/chem_environment_config.h:25` (`42 mmHg ~ 5.5e-2 mol/m^3`), `1 mmHg (1 Torr) = 1.3095e-3 mol/m^3`, `1 uM = 1.0e-3 mol/m^3`, and `% O2` at 1 atm dry `= pct * 7.60 * 1.3095e-3`.

| ID | Paper (Year, DOI) | Depth | System / Assay | Raw Measurement | GutIBM Value | Discrepancy vs. Baseline |
| --- | --- | --- | --- | --- | --- | --- |
| O2-01 | Špela Konjar et al. (2021, `10.3390/ijms22179170`) | Full-Text | Mammalian GI longitudinal mucosa (review) | 58 -> 32 -> 11 -> 3 mmHg | 0.0144 asc; 3.93e-03 sig (mol/m^3) | Shipped 55.0e-6 is ~71x below sigmoid (3 mmHg) and ~262x below ascending colon (11 mmHg) |
| O2-02 | Alexander V. Zhdanov et al. (2016, `10.1007/s00018-016-2323-x`) | Abstract-Only | Mouse distal colon ex vivo (top 50-60 um mucosa) | 5-35 (4-28) uM (Torr) | 5.00e-03 to 0.035 (mol/m^3) | Shipped 55.0e-6 is ~91x to 636x below normal mucosal 5-35 uM (0.005-0.035 mol/m^3) |
| O2-03 | Lindsey G. Albenberg et al. (2014, `10.1053/j.gastro.2014.07.020`) | Blocked-In-FullText | Mouse cecum lumen & mucosal interface | <1 mmHg | <1.31e-03 (mol/m^3) | Confirms lumen floor <1.31e-3 mol/m^3 while host vasculature/mucosa drives radial gradient |
| O2-04 | Raehyun Kim et al. (2019, `10.1088/1758-5090/ab446e`) | Full-Text-Chunks | Primary human colonic epithelium + anaerobes | 0.1-1 mmHg (<1% O2); 0.8+/-0.1% vs 11.1+/-0.5% % O2 / mmHg | 1.31e-04-1.31e-03 phys lum; 7.96e-03 dev lum; 0.11 basal (mol/m^3) | Basal epithelial supply ~1.10e-1 mol/m^3 (84 mmHg) vs physiological lumen 1.31e-4 to 1.31e-3 mol/m^3 (device lumen 7.96e-3 mol/m^3) |
| O2-05 | Jiaquan Yu et al. (2025, `10.64898/2025.12.16.694730`) | Blocked-In-FullText | Human/mammalian intestinal mucosal barrier | 10-60 uM/mm | 0.01 to 0.06 per mm (mol/m^3/mm) | Over 200-500 um mucus layer, supports ~2e-3 to 3e-2 mol/m^3 drop from epithelium to lumen |

### 3.2 Anaerobic Maintenance Factor (`chem_env.oxygen.anaerobic_maintenance_factor`)

In `src/fixes/fix_metabolism.cpp:662-667`, `anaerobic_maintenance_factor` scales `cfg_.carbon_maintenance_rate` linearly with `realized_fermentation_fraction`. Biothermodynamically, this equals `(m_ATP,anaer / m_ATP,aer) * (Y_ATP/Glc,aer / Y_ATP/Glc,anaer)`.

| ID | Paper (Year, DOI) | Depth | System / Assay | Raw Measurement | GutIBM Value | Discrepancy vs. Baseline |
| --- | --- | --- | --- | --- | --- | --- |
| MAINT-01 | Karel Olavarría et al. (2019, `10.1007/s00253-019-09909-6`) | Full-Text | E. coli K-12 MG1655 (aer) & BW25113 (anaer) | 0.054 anaer vs 0.075+/-0.015 aer mol ATP / C-mol_BM / h | 0.72 (ATP-basis ratio) (dimensionless) | On an ATP basis, anaerobic m_ATP is slightly LOWER (0.72x), NOT higher; high carbon cost comes purely from lower ATP/Glc yield |
| MAINT-02 | Karel Olavarría et al. (2019, `10.1007/s00253-019-09909-6`) | Full-Text | E. coli K-12 MG1655 (aer, D=0.049) vs GJT001 (anaer, D=0.1) | 3.55 anaer vs 0.77 aer mmol Glc / gCDW / h | 4.61 (dimensionless ratio) | 4.61x total q_Glc ratio (with Y_BM 0.156 vs 0.426 g/g) is ~3.3x below shipped 15.0 |
| MAINT-03 | Simon Boecker et al. (2021, `10.15252/msb.202110504`) | Full-Text | E. coli K-12 MG1655 anaerobic glucose WT (n=3) | 4.41+/-0.54 (growing); 4.16+/-0.27 (arrested) mmol ATP / gDW / h | 4.16 - 4.41 (controls 3.30-5.42) (mmol ATP / gDW / h) | Supports ~2.6-2.9 ATP/Glc mixed-acid fermentation and moderate NGAM (~4.2-4.4 mmol ATP/gDW/h), inconsistent with 15x carbon maintenance |
| MAINT-04 | Ross P. Carlson et al. (2026, `10.1002/1873-3468.70380`) | Full-Text | E. coli aerobic respiratory chain | 1.0-1.5 (exp 1.01-1.07 at mu_max) mol ATP / O atom | 3.36 to 5.76 (w/ 0.72 m_ATP); 4.67-8.00 (ATP yield) (dimensionless factor) | Shipped 15.0 is ~2.6-4.4x too high; realistic E. coli carbon maintenance factor is ~3.4-5.8 (upper bound ~8.0 if m_ATP equal) |
| MAINT-05 | Walter P. Hempfling & Stanley E. Mainzer (1975, `10.1128/jb.123.3.1076-1087.1975`) | Blocked-In-FullText | E. coli B (glucose-limited chemostat) | 0.9 mg-atoms O/gDW/h (aer); 18.9 mmol ATP/gDW/h (anaer) mg-atoms O or mmol ATP/gDW/h | P/O=1.5..2 -> m_ATP_aer=1.35..1.8; strain B m_ATP ~4x higher than K-12 (mmol ATP / gDW / h) | E. coli B exhibits much higher anaerobic m_ATP (18.9) than K-12 MG1655 (4.16-4.41); not representative of K-12 baseline |
| MAINT-06 | A Varma & Bernhard Ørn Palsson (1994, `10.1128/aem.60.10.3724-3731.1994`) | Blocked-In-FullText | E. coli W3110 | 10.5 aer vs 18.5 anaer; NGAM=7.6 mmol/gDW/h | Max q_Glc ratio = 1.76x; NGAM=7.6 mmol ATP/gDW/h (dimensionless / mmol ATP/gDW/h) | Anaerobic q_Glc capacity is 1.76x aerobic; supports ~4-6x maintenance carbon factor (comparable to anaerobic_carbon_cost_factor=4.1) |

**Calibration conclusion:** Set `chem_env.oxygen.epithelial_conc` in the colonic range **`5.0e-3 to 1.5e-2 mol/m^3`** (`~4–11 mmHg`, treating `42 mmHg = 5.5e-2 mol/m^3` as an arterial/proximal upper bound and `55.0e-6` as a confirmed 1,000x unit bug), and lower `chem_env.oxygen.anaerobic_maintenance_factor` from `15.0` to **`4.5` (sweep `3.4–5.8`)**.

---

## 4. Pillar 3: Mammalian GI Transit Times and Colicinogeny Turnover

GutIBM defaults to `layer3.transit_h = 12.0 / 8.0 / 10.0 h` (`30.0 h` total).

| ID | Paper (Year, DOI) | Depth | System / Assay | Raw Measurement | GutIBM Value | Discrepancy vs. Baseline |
| --- | --- | --- | --- | --- | --- | --- |
| TR-01 | Nicola Procházková et al. (2024, `10.1038/s41564-024-01856-x`) | Full-Text | Human healthy adults (n=50) | CTT 13.6 (2.1-63.5); WGTT 23.3 (12.4-72.3) h | 13.6 CTT (range 2.1-63.5); 23.3 WGTT (h) | Shipped 30.0 h sits inside observed range (2.1-63.5 h) but is 2.2x slower than SmartPill healthy median (13.6 h) |
| TR-02 | Gursharan Kaur Nandhra et al. (2023, `10.3390/jcm12165272`) | Full-Text | Human healthy adults (22 ingestible capsule studies) | CTT 15.9-28.9; WGTT 23.0-37.4 h | 15.9-28.9 CTT; 23.0-37.4 WGTT (h) | Shipped 30.0 h matches upper boundary of IC study medians (15.9-28.9 h) and lower-mid ROM colonic norms (~30-70 h) |
| TR-03 | Gursharan Kaur Nandhra et al. (2019, `10.1111/nmo.13734`) | Blocked-In-FullText | Human healthy adults (n=111) | F: 32/34/17/17%; M: 33/25/14/28% % of CTT | At 20h CTT: ~6.5/5.9/3.1/4.5 h (Asc/Trans/Desc/Rectosig) (h (scaled)) | Confirms ascending/cecal colon holds the largest share (~32-33%) among proximal/mid/distal segments |
| TR-04 | Chloë Rotsaert et al. (2025, `10.1093/jas/skaf107`) | Abstract-Only | Piglets (n=10, TopigsNorsvin x Pietrain) | 12.75+/-1.7 (12h45m+/-1h42m) h | 12.75+/-1.7 (h) | Young piglet total transit (~12.75 h) is 2.35x faster than human colonic default (30 h) |
| TR-05 | Laura J. Henze et al. (2020, `10.1016/j.ejps.2020.105627`) | Full-Text-Chunks | Landrace pigs (fed/fasted) | SITT 2.3-4.0; GET 20-233 (up to 264) h | SITT 2.3-4.0 (gastric retention 20-264) (h) | Rigid large capsules show massive gastric delay (20-264 h) in pigs; SITT (2.3-4.0 h) tracks human SI |
| TR-06 | Andre Tan et al. (2016, `10.1016/j.livsci.2016.03.010`) | Abstract-Only | Young pigs (piglets) | Dye: 13.2 first, 24.1 last; ROM: 25.2 & 48.9 h | 13.2-24.1 fluid; 25.2-48.9 solid (h) | Pig solid WGTT (25-49 h) overlaps adult human WGTT; fluid phase (13-24 h) is ~2x faster |
| TR-07 | Andre Tan et al. (2018, `10.1111/ner.12761`) | Abstract-Only | Young pigs | WGTT 7.7-72.2 (med 23.0-28.4); bowel+rectum 5-53 h | 5-53 bowel+rectum; 23.0-28.4 WGTT med (h) | Large-bowel transit span (5-53 h) brackets the 30.0 h human colonic parameterization |
| TR-08 | Eva A. Latymer et al. (1990, `10.1080/17450399009428415`) | Abstract-Only | Growing/adult pigs (30-85 kg) | Colon: 4.5-50.3 (peak 16.5-48.8); rectum: 24-<51 h | Colon peak 16.5-48.8; rectum 24-<51 (h) | Growing/adult pigs closely match human colonic transit timescales (~20-48 h) |
| TR-09 | Gaku Wagai et al. (2024, `10.1007/s10620-024-08420-x`) | Full-Text | Inbred mice (C57BL/6, BALB/c, C3H, DBA/2J) | WGTT 2-6 (C57BL/6, BALB/c); distal colon <3 (C3H); MRT 10.3 (DBA/2J) h | 2-6 WGTT (colonic ~1-3 h) (h) | Mouse colon transit (~1-3 h) is ~10x faster than GutIBM human default (30 h); must rescale when modeling mouse experiments |
| TR-10 | Parasuraman Padmanabhan et al. (2013, `10.1186/2191-219x-3-60`) | Full-Text | Mice (in vivo SPECT imaging) | 1h stomach+SI; 3h cecum; 6-8h colon/rectum (~6h WGTT) h | ~3-5 colonic residence (6h WGTT) (h) | Supports 3-5 h cecal+colonic transit for mouse simulations (6-10x faster than human 30 h) |
| TR-11 | Sean T. Koester et al. (2021, `10.1016/j.xpro.2021.100938`) | Full-Text | Conventional & germ-free mice | 2-8 (conventional shorter; GF longer) h | 2-8 WGTT (h) | Independently confirms 2-8 h murine WGTT and shows germ-free status slows transit |
| TR-12 | Belinda Barnes et al. (2007, `10.1099/mic.0.2007/007120-0`) | Blocked-In-FullText | Australian mammals (carnivores vs herbivores/omnivores) | Carnivore < herbivore/omnivore colicin freq. qualitative / model | High turnover (short transit) selects against colicinogeny (ecological scaling law) | Directly links Pillar 3 transit time to Pillar 4 ColE2 selection: fast mouse transit (~2-6 h) suppresses producer invasion unless crypt/mucus retention holds |

**Calibration conclusion:** Keep `30.0 h` (`12/8/10 h`) as the standard adult human (and growing-pig) baseline while adding a fast-transit human sensitivity arm (`14–20 h` matching SmartPill median `13.6 h`), and introduce a **`1.0 / 1.5 / 1.5 h` (`4.0 h` colonic; `2–6 h` WGTT)** preset whenever validating against murine colonization data.

---

## 5. Pillar 4: ColE2 Exclusion Dynamics & Prophage Induction Priors

### 5.1 ColE2 Invasion, Lysis Gating, and Spatial Benchmarks

| ID | Paper (Year, DOI) | Depth | System / Assay | Raw Measurement | GutIBM Value | Discrepancy vs. Baseline |
| --- | --- | --- | --- | --- | --- | --- |
| COL-01 | Erik Bakkeren et al. (2025, `10.1038/s41564-025-02162-w`) | Full-Text | Mouse gut (Salmonella / E. coli + 4% sorbitol) | 0% exclusion by toxin alone; 72h exclusion w/ private nutrient hours / qualitative | 72h exclusion only when private niche + ColE2 act together (ecological rule) | Explains why pure ColE2 cannot invade high-B12/high-overlap settings; requires low resource overlap + active colicin |
| COL-02 | Felix J.H. Hol et al. (2014, `10.1186/s12915-014-0068-2`) | Full-Text | E. coli pColE2 microfluidic & batch cultures | M9-Glc: 22-50% Psos, 1% lyse (n=548); LB: 1.6-26% Psos, 44%/15% lyse (n=815) % induced; % lysed | 0.01 (M9-Glc lysed/induced) vs 0.15-0.44 (LB) (probability) | Shows lysis probability is nutrient-regulated via CsrA (1% in M9-Glc vs 15-44% in LB); relieving CsrA (pColE2-TT) collapses diversity |
| COL-03 | Anna S. Weiß et al. (2020, `10.1038/s41598-020-61086-z`) | Full-Text | E. coli ColE2 (CWT, CAMP, CREP1/2) +/- MitC | 70-85 min (CREP1/2) -> 120-150 min (CAMP/CWT); R^2=0.89 min; R^2 | 4200-9000 s post-SOS lysis delay; S wins 1:100 at 0 MitC (s; outcome) | At 1:100 C:S without external DNA stress, S excludes C; longer cea-cel delay (19 min vs 6 min; 150 min lysis) boosts toxin yield |
| COL-04 | Lusine Ghazaryan et al. (2019, `10.3390/microorganisms7110564`) | Full-Text | E. coli BZB1011 ColE2 vs WT-GFP / ColA (24h) | w=1.94+/-0.38 (ColE2); threshold >0.5-2% (vs WT), >1-10% (vs ColA) w (dimensionless); % init | ColE2 4x basal expression vs ColE7; threshold 0.005-0.02 (proportion) | Benchmark for Spec 14 V1: ColE2 invades WT at >0.5-2% initial fraction with w=1.94+/-0.38 over 24h when free B12 is not saturating |
| COL-05 | Thierry Kuhn et al. (2023, `10.1098/rstb.2021.0503`) | Full-Text | E. coli producer/filtrate + lattice RPS PDE | >40% S CFU drop in 30 min; D=0.002, dl=3e-5..3e-4, u=0.0075..0.03 % kill / 30 min; dimensionless | k_kill > 2.84e-4 /s at 50% filtrate; intermediate retention required (1/s; rule) | Validates intermediate toxin retention regime: washout or excessive immobilization both break ColE2 killing |
| COL-06 | Benjamin Kerr et al. (2002, `10.1038/nature00823`) | Abstract-Only | E. coli C-S-R triad on agar plates / flasks | Coexistence on static plate; extinction in mixed flask qualitative benchmark | Local dispersal preserves C-S-R coexistence (validation target) | Core in vitro spatial validation benchmark: local dispersal maintains all three strains whereas global mixing collapses diversity |
| COL-07 | Benjamin C Kirkup & Margaret A. Riley (2004, `10.1038/nature02429`) | Blocked-In-FullText | Mouse gut (streptomycin-treated C-S-R triad) | Coexistence of colicin producer, sensitive & resistant in mouse gut qualitative benchmark | Local gut structure preserves C-S-R coexistence (validation target) | Core in vivo spatial validation benchmark: mouse gut spatial structure sustains RPS coexistence |

### 5.2 Prophage Induction Priors and Blocked-Source Transparency

| ID | Paper (Year, DOI) | Depth | System / Assay | Raw Measurement | GutIBM Value | Discrepancy vs. Baseline |
| --- | --- | --- | --- | --- | --- | --- |
| PRO-01 | Caroline Henrot et al. (2026, `10.64898/2026.05.26.727859`) | Blocked-In-FullText | 1,768 human fecal E. coli + dixenic mice | 6h in vitro exclusion -> 10d mouse coexistence via lysogens h; days | 6h (21600 s) in vitro vs 10d (864000 s) lysogen coexistence (s) | Only DNA-damaging bacteriocins (E-type endonuclease colicins like ColE2 + Mcc1229) trigger prophage induction and lysogen conversion |
| PRO-02 | Marianne De Paepe et al. (2016, `10.1371/journal.pgen.1005861`) | Full-Text | Monoxenic germ-free mice (E. coli lambda, 6 mice) | 0.016 h^-1 (ODE); 0.0170+/-0.0046 /gen (reg) vs 3.4e-4 /gen in vitro h^-1; per gen | 4.44e-6 /s (ODE h^-1); ~1.18e-6 /s at d=0.25 h^-1 (~50x in vitro) (1/s) | De Paepe's raw mouse rate (~1.18e-6 to 4.44e-6 /s) is 7-25x HIGHER than Spec 14's 1.75e-7 /s (1.5%/day), which conflates /gen with /day |
| PRO-03 | Nelson Frazão et al. (2022, `10.1038/s41467-022-33412-8`) | Full-Text | Mouse gut evolved E. coli lysogens (d165-d493) | 0.67-0.69 h^-1 (new) -> 0.37-0.39 h^-1 (evolved); 0.41 (resident) h^-1 | 1.03e-04 to 1.92e-04 (1/s) | Provides empirical upper bound on SOS-induced lysis rate (~1.03e-4 to 1.92e-4 /s) and shows gut adaptation halves inducibility |
| PRO-04 | Jie Hu et al. (2023, `10.1038/s41522-023-00445-w`) | Full-Text | Mice + E. coli ATCC 25922 (Phi1-Phi3) + Clostridia | 2.13x in vitro; 11x in vivo (3.59e9 -> 3.98e10 copies/100 mg) fold change; copies/100 mg | 2.13x in vitro; 11x in vivo (dimensionless multiplier) | Shows physiological microbiota metabolites (propionate 2.5-7.5 mM) modulate RecA-dependent prophage induction up to 11x in vivo |
| PRO-05 | Yan Huang et al. (2025, `10.1016/j.cub.2025.01.006`) | Abstract-Only | Human gut microbiome (corroborating 10.1016/j.cub.2025.03.073 [Blocked]) | 0.001 - 0.01 events / bacterium / day | 1.16e-08 to 1.16e-07 (realized) (1/s) | Spec 14 default (0.0120-0.0121/bact/day at prev=0.8; 0.0151 at 1.0) EXCEEDS 0.01/bact/day ceiling before adding Hill damage induction |

**Transparency note on `10.1016/j.cub.2025.03.073` and upstream metadata:** Full-text retrieval for Lopez et al. (*Current Biology*, `10.1016/j.cub.2025.03.073`) returned `status = "blocked"` from the publisher. The quantitative human gut community induction estimate (`0.001–0.01 events per bacterium per day`, corresponding to `-ln(1 - p)/86400 = 1.16e-8 to 1.16e-7 /s` realized) is verified verbatim in the accompanying Dispatch abstract by Huang, Wang & Ma (2025, *Current Biology*, `10.1016/j.cub.2025.01.006`). Additionally, note that `10.1038/nature02429` (Kirkup & Riley 2004, *Nature* vol. 428) carries an upstream Consensus `journal_name` tag of `Nature Cell Biology`; the correct journal (`Nature`) is noted alongside the upstream tag in `consensus_empirical_benchmarks.csv`.

---

## 6. Recommended Cycle 2 Parameter Updates

| Parameter Key | Current / Proposed | Recommended Cycle 2 Value / Sweep | Primary Empirical Grounding |
| --- | --- | --- | --- |
| `b12.initial_conc` (free/in vitro) | `1.0e-3 mol/m^3` | **`1.0e-6 to 1.5e-5 mol/m^3`** (`3.69e-6` MacFarlane; `1.35e-5` human true cobalamin) | Allen & Stabler 2008; Kelly 2019; Forgie 2023; Kundra 2022 |
| `chem_env.oxygen.epithelial_conc` | `55.0e-6 mol/m^3` | **`5.0e-3 to 1.5e-2 mol/m^3`** (`4–11 mmHg`; `5.5e-2` upper bound) | Konjar 2021; Zhdanov 2016; Albenberg 2014; Kim 2019 |
| `chem_env.oxygen.anaerobic_maintenance_factor` | `15.0` | **`4.5`** (sweep **`3.4–5.8`**) | Olavarría 2019; Boecker 2021; Carlson 2026 |
| `layer3.transit_h` (human vs. mouse) | `12/8/10 h` (`30.0 h`) | Human: `12/8/10 h` + fast arm `6/4/4 h`; **Mouse: `1.0/1.5/1.5 h` (`4.0 h`)** | Procházková 2024; Nandhra 2023; Wagai 2024; Padmanabhan 2013 |
| `spec14.spontaneous_induction_rate` | `1.75e-7 /s` at `prev=0.8` | Constrain `P(inducible) * rate` in **`1.16e-8 to 1.16e-7 /s`** (`0.001–0.01 /d`); test `~1e-6 /s` as monoxenic mouse bound | Huang 2025 (`10.1016/j.cub.2025.01.006`); De Paepe 2016; Frazão 2022 |
