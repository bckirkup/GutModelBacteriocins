Question: Provide literature-grounded parameter estimates for modeling E. coli 
    colicin-mediated rock-paper-scissors (Producer-Sensitive-Resistant) dynamics 
    in the colonic mucus layer. This is for an individual-based model (IBM) with 
    explicit spatial structure, diffusion, and receptor-mediated colicin killing.

    I need quantitative estimates with citations for ALL of the following:

    ## 1. RESISTANCE FITNESS COST

    What is the growth rate cost of colicin resistance in E. coli?

    - Resistance typically arises from loss/modification of outer membrane 
      receptors (BtuB for colicin E, FepA for colicin B, CirA for colicin Ia)
    - These receptors have primary functions: BtuB = vitamin B12 uptake, 
      FepA = ferric enterobactin uptake, CirA = ferric catecholate uptake
    - Loss of these receptors = loss of nutrient uptake = growth cost
    - Need: percent reduction in max growth rate (mu_max) for BtuB-null, 
      FepA-null, and double-null mutants
    - Context: colon environment with available iron and B12
    - Key authors: Feldgarden, Riley, Kerr, Chao, Levin, Gordon
    - Look for competition experiments measuring relative fitness w/w0

    ## 2. COLICIN IMMUNITY PROTEIN EXPRESSION COST

    Producers carry the immunity gene. What is the metabolic cost of 
    constitutive immunity protein expression?
    - Typically small (immunity protein is small, ~10-15 kDa)
    - But it is expressed constitutively from its own promoter
    - Need: any measured fitness cost of carrying the immunity gene alone
      (without the colicin structural gene)

    ## 3. COLONIC CORRINOID CONCENTRATIONS

    What is the concentration of corrinoids (vitamin B12 analogs) in the 
    human colonic mucus layer?
    - Total corrinoids vs true cobalamin (cyanocobalamin/methylcobalamin)
    - Important because non-cobalamin analogs may compete for BtuB binding 
      but with different affinity
    - Need: mol/L or mol/m3 concentrations
    - Sources: human fecal metabolomics, colonic fluid analysis, 
      gut microbiome B12 production studies
    - Key distinction: luminal vs mucosal concentrations

    ## 4. BtuB-CORRINOID BINDING AFFINITY (Kd)

    What is the Kd for:
    - Cobalamin binding to BtuB: literature value
    - Colicin E binding to BtuB: literature value  
    - Non-cobalamin corrinoid analogs binding to BtuB: if available
    - Competitive inhibition: does B12 compete with colicin for BtuB?
    - Key papers: Cadieux, Bhatt, Kurisu, Cherezov (structural/binding studies)

    ## 5. IMMIGRATION / COLONIZATION RATE

    What is the rate of new E. coli strain arrival in the human colon?
    - From diet, environment, other hosts
    - Strain turnover rate in longitudinal studies
    - Need: events per day or per week
    - Key studies: Caugant, Poulsen, Leatham, Conway, Gordon longitudinal studies
    - Relevant: how many distinct E. coli strains coexist in a healthy gut?

    ## 6. E. COLI POPULATION DENSITY IN COLONIC MUCUS

    What is the in vivo density of E. coli in the mucus-associated community?
    - Need: CFU per gram mucus or per cm2 mucosal surface
    - Distinguish: mucus-associated vs luminal/fecal
    - Healthy vs dysbiotic ranges
    - Spatial patchiness: colony sizes, inter-colony distances

    ## 7. COLICIN DIFFUSION IN MUCUS

    What is the diffusion coefficient of colicin proteins in colonic mucus?
    - Colicins are ~40-70 kDa proteins
    - Mucus gel retards diffusion relative to free solution
    - Need: D_free in water AND retardation factor in mucus
    - Any measurements of protein diffusion in reconstituted or native mucus

    ## 8. SOS INDUCTION AND LYSIS PARAMETERS

    - Spontaneous SOS induction rate in E. coli (per cell per generation)
    - Fraction of SOS-induced cells that lyse (vs repair and survive)
    - Burst size: number of colicin molecules released per lysis event
    - Time from induction to lysis
    - Key papers: Cascales, Lazdunski, Riley, Gillor

    ## 9. COLICIN KILLING KINETICS

    - Is colicin killing truly "one-hit" (single molecule sufficient)?
    - What is the minimum number of molecules needed to kill a cell?
    - How fast is killing once the colicin binds the receptor?
    - Temperature and pH dependence relevant to the colon (37C, pH 6-7)

    ## 10. CONJUGATION RATE FOR COLICIN PLASMIDS

    - Transfer rate for Col plasmids (ColE1, ColB) in structured environments
    - ColE1 is typically non-conjugative (mobilizable only)
    - ColB plasmids may be conjugative
    - Need: transfer events per donor-recipient pair per hour
    - In biofilm/mucus vs liquid culture

    For each parameter, provide:
    - Best estimate value with units
    - Range (min-max from different studies)
    - Key citation(s) with full BibTeX
    - Caveats about applicability to the colonic mucus environment

    Target: 30-50 citations.

# Parameterization of colicin-mediated Producer–Sensitive–Resistant dynamics in colonic mucus

## Evidence grading

**D** denotes a direct measurement of the requested quantity; **A**, a close analogue; **I**, an inferred or weakly informed modeling prior. Several requested parameters—especially receptor-specific growth costs in colon-like medium, immunity-only cost, mucosal corrinoids, colicin burst size, and ColB conjugation rate—have not been measured directly. Those values should therefore be varied in sensitivity analyses rather than treated as established constants.

## Recommended IBM parameter set

| Process | Best estimate | Sensitivity range | Evidence and caveat |
|---|---:|---:|---|
| **BtuB-null growth cost** | **2% reduction in μmax** | 0–10% | **I.** BtuB loss is often effectively neutral when B12-dependent metabolism is dispensable. Apply a larger cost only under corrinoid-dependent methionine synthesis or ethanolamine/1,2-propanediol utilization. Feldgarden and Riley demonstrated that colicin-resistance phenotypes can carry fitness effects, but the available literature does not establish a universal BtuB-null percentage. Use `μR=0.98 μS` as a neutral-centered prior, not a measured colon value. |
| **FepA-null growth cost** | **5%** | 0–20% | **I.** Iron uptake is redundant; FepA loss may be nearly neutral in iron-replete medium but costly under iron limitation. A single-receptor deletion should therefore have a smaller cost than loss of multiple catecholate receptors. |
| **BtuB/FepA double-null cost** | **8%** | 0–30% | **I.** Use approximately multiplicative costs, `1-(1-cBtuB)(1-cFepA)`, unless iron/corrinoid limitation is modeled explicitly. The upper tail represents nutritional immunity or strong dependence on the lost pathway. |
| **CirA-null cost** | **3%** | 0–15% | **I.** CirA and Fiu/FepA functions overlap. In iron-replete conditions use near zero; under catecholate-limited conditions increase it. |
| **Immunity protein alone** | **0.2% μmax cost** | 0–2% | **I.** No convincing immunity-only competition measurement was found. Small 10–15-kDa proteins expressed at low constitutive levels should have a cost below ordinary assay resolution. A ColIb plasmid comparison found no measurable growth-rate effect, although that does not isolate immunity expression. |
| **Total fecal corrinoids** | **0.97 μM** | 0.2–2.4 μM | **D→derived.** Mean 1309 ng/g wet feces divided by ≈1350 g/mol and assuming 1 kg/L. Human mucus was not measured. |
| **True cobalamin** | **14 nM** | 3–100 nM | **D→derived.** True Cbl was 19 ng/g, 1.4% of total; ≈98% consisted of analogues. |
| **Non-Cbl analogues** | **0.95 μM** | 0.2–2.3 μM | **D→derived.** Dominant species was 2-methyladeninyl cobamide, 60.6%; p-cresolyl 16.3%, adenyl 12.5%, and cobinamide 1.8%. |
| **Mucus corrinoid concentration** | **0.3 μM total; 10 nM Cbl** | 0.03–1 μM total | **I.** No human colonic-mucus measurement was located. Use a mucus:lumen partition coefficient of 0.3 (0.03–1) and calibrate independently. |
| **BtuB–cobalamin Kd** | **0.3 nM** | 0.05–5 nM | **Literature-based I.** High-affinity TonB-dependent uptake is well established, but published values vary with membranes, detergent, Ca²⁺ and assay. Use receptor occupancy rather than assuming every fecal corrinoid is freely available. |
| **BtuB–colicin E receptor-domain Kd** | **1 nM** | 0.1–20 nM | **Literature-based I.** Structural and thermodynamic studies support nanomolar binding. Values are colicin- and construct-dependent; full-length toxin and isolated receptor domains need not agree. |
| **Analogue–BtuB Kd** | **10 nM** | 0.3 nM–1 μM | **I.** Analogue-specific BtuB constants are sparse. Assign each analogue a relative affinity and fit it; do not treat total corrinoids as cobalamin equivalents. |
| **B12 competition with colicin E** | **competitive receptor occupancy** | explicit 0–100% displacement | Vitamin B12 can protect/rescue cells during early BtuB-dependent colicin binding, supporting competition or steric interference. Implement `fCol = ([Col]/KdCol)/(1+[Col]/KdCol+Σ[Ci]/Kd,i)`; subsequent OmpF/Tol or TonB translocation remains separate. |
| **Candidate strain exposure** | **5 strains host⁻¹ day⁻¹** | 1–20 day⁻¹ | A recent synthesis estimates 1–20 strains/day may contend for colonization. This is exposure, not successful establishment. |
| **Successful immigration** | **0.005 day⁻¹ per modeled patch** | 10⁻⁴–0.1 day⁻¹ | **I.** Multiply exposure by a small establishment probability and scale by simulated surface area. Do not inject 1–20 established strains/day into every microscopic patch. |
| **Concurrent dominant strains** | **2** | 1–4 | Historical studies generally recovered one or two serotypes per time point, while at least four could overlap. Limited colony sampling undercounts rare strains. |
| **Transient/resident persistence** | **7 d / 365 d** | 3–30 d / 30–1065 d | Rapid and slow turnover occur on approximately week and year scales; a 20-y series identified 35 clones and minimum carriage durations of 7–1065 d. |
| **Mucus-associated E. coli density** | **10⁷ cells g⁻¹ mucus** | 10⁴–10⁹ g⁻¹ | **I.** A direct healthy-human mucus-specific value is unavailable. Fecal/Enterobacteriaceae abundance spans 10²–10⁹ CFU/g. Use 10⁷ only as a central colonization-density prior; dysbiosis can approach the upper bound. |
| **Equivalent areal density** | **10⁵ cells cm⁻²** | 10²–10⁷ cm⁻² | **Derived.** Assumes a 100-μm layer over 1 cm², density 1 g/cm³, and 10⁷ cells/g. Patchiness should be generated explicitly rather than imposed as uniform density. |
| **Initial microcolony radius** | **10 μm** | 2–100 μm | **I.** Direct healthy-colon inter-colony distances for E. coli were not found. Use clustered seeding and test center-to-center spacings of 20–500 μm. |
| **Colicin D in water, 25°C** | **8.5×10⁻¹¹ m² s⁻¹** | 7–10×10⁻¹¹ | **A.** Ovalbumin, F(ab)/Fc, and BSA (45–68 kDa) had 8.3–8.7×10⁻¹¹ m²/s in water. |
| **Colicin D in mucus, 37°C** | **7×10⁻¹¹ m² s⁻¹** | 2×10⁻¹¹–1.2×10⁻¹⁰ | **A/I.** Size-matched proteins in fresh cervical mucus at 25°C had 5.7–8.8×10⁻¹¹ m²/s and `Dmucus/Dwater=0.68–1.0`. Colonic mucus may bind colicins more strongly; include reversible mucin binding or a low-D subpopulation. |
| **Spontaneous colicin/SOS-on fraction** | **1% per generation** | 0.5–3% | **A.** Spontaneous reporter-positive fractions across colicins were ≈0.5–3%; a ColIb multicopy reporter gave 0.93%. This is an ON fraction, not necessarily a transition probability; for a one-generation ON lifetime they are approximately equal. |
| **Fraction of induced producers that lyse** | **0.2** | 0.05–1 | Under strong DTPA+mitomycin-C induction, ≈80% expressed ColIb but 20.13% expressed prophage lysis genes. Group-A colicins with dedicated lysis genes may differ markedly. |
| **Induction-to-lysis delay** | **90 min** | 30–240 min | **I.** Strongly system-, growth-, and stress-dependent. Use a gamma/lognormal delay, not an exponential event if synchronized bursts matter. In one induced ColIb/prophage system population decline occurred after ≈2 h. |
| **Colicin burst per lysed cell** | **10⁵ molecules** | 10³–10⁶ | **I; major data gap.** No direct single-cell molecule count was found. This range corresponds to roughly 0.07–70 fg for a 45-kDa toxin and should be calibrated against inhibition-zone or liquid-killing data. |
| **Lethal-hit threshold** | **1 successfully translocated toxin** | 1–10 | Classical “one-hit” kinetics support one productive toxin event, not necessarily one extracellular molecule. Receptor binding, failed translocation, sequestration, and immunity make the extracellular molecules-per-death much larger. |
| **Binding-to-commitment time** | **5 min** | 1–20 min | **I/A.** Classical adsorption experiments quantify substantial killing within ≈10 min at 37°C; one study explicitly incubated toxin-cell mixtures for 10 min (retrieved de Zwaig and Luria paper). Use separate reversible binding and irreversible translocation steps. |
| **Commitment-to-loss-of-viability** | **15 min** | 2–60 min | Pore formers can depolarize rapidly; nuclease colicins require translocation and catalytic damage. Population-level ColIb effectiveness fitted at 4.95 OD⁻¹ h⁻¹ at 37°C provides a useful calibration target. |
| **pColE1 self-transfer** | **0** | exactly 0 without helper | ColE1 lacks a complete conjugation apparatus. Model mobilization only if a compatible conjugative helper supplies transfer functions. |
| **pColE1 mobilization with helper** | **10⁻¹¹ mL cell⁻¹ h⁻¹** | 10⁻¹³–10⁻⁹ | **I.** Highly helper-, host-, and contact-dependent. Prefer a contact-event probability in the IBM rather than a well-mixed mass-action rate. |
| **Conjugative ColB/ColV-like plasmid** | **10⁻¹⁰ mL cell⁻¹ h⁻¹** | 10⁻¹²–10⁻⁹ | **I.** No defensible direct ColB pairwise rate was found. Biofilms can increase transfer by orders of magnitude in some plasmid-host combinations, but architecture may also restrict donor-recipient interfaces. Use a structured-environment multiplier of 0.1–100, fitted to the chosen plasmid. |

## Interpretation and implementation notes

### Resistance costs

There is no literature basis for assigning a single universal cost to “colicin resistance.” Receptor loss is condition dependent: BtuB matters only when a B12-dependent pathway limits growth, whereas FepA/CirA costs increase under iron limitation and can be masked by redundant siderophore receptors. The colon is not simply “iron available”: host nutritional immunity and microbial siderophore competition can make freely accessible ferric iron scarce. Implement costs mechanistically where possible, for example

`μ = μbase × min(gC, gN, gFe, gCorrinoid)`,

with receptor genotype changing only the corresponding uptake term. If a fixed-cost RPS model is required, the central estimates above preserve the expected ordering `BtuB < FepA < double-null` but must be subjected to broad sensitivity analysis. Feldgarden and Riley’s resistance studies are the key historical evidence, but they do not provide colon-specific isogenic μmax penalties.

### Corrinoid conversion

The fecal mean corresponds to approximately `1.309 mg L⁻¹ / 1350 g mol⁻¹ = 0.97 μmol L⁻¹ = 0.97 mmol m⁻³`. True Cbl is approximately `0.019 mg L⁻¹ / 1355 g mol⁻¹ = 14 nmol L⁻¹ = 0.014 mmol m⁻³`. These are wet-feces concentrations, not free aqueous activities. Binding to cells/solids, intracellular corrinoids, and mucus partitioning can reduce the freely competing concentration substantially. Allen and Stabler found the dominant material to be noncanonical cobamides rather than true Cbl.

### Spatial scaling

At `D=7×10⁻¹¹ m²/s`, the root-mean-square displacement is approximately 0.71 mm in one hour in three dimensions, before binding or degradation. Thus spatial localization in an IBM will arise primarily from finite lifetime, adsorption to receptors/debris/mucin, and producer/sensitive geometry—not from gel obstruction alone. The cervical-mucus analogue was ≈95% water and allowed many proteins to diffuse at 0.7–1.0 of their water value. Native inner colonic mucus can differ substantially in mucin concentration and chemistry.

### Killing kernel

A mechanistic sequence is preferable:

1. reversible receptor binding with nutrient/corrinoid competition;
2. TonB- or Tol-dependent translocation;
3. irreversible lethal commitment with probability `ptrans`;
4. a mode-specific death delay.

“One hit” should be applied at step 3. It does not imply that every molecule released kills one bacterium. O-antigen and growth conditions can occlude outer-membrane receptors, so receptor number alone is insufficient.

### Immigration

The 1–20 strains/day estimate represents exposure to the whole host, whereas successful colonization is much rarer. A useful IBM decomposition is

`λpatch = λexposure × pestablish × Apatch/Acolon`.

Longitudinal data support both short-lived transients and persistent residents. In one 20-y series, 210 isolates yielded 35 clones: 25 subdominant/transient, five dominant/transient, and five dominant/resident. Sampling only 10–15 colonies per stool misses low-frequency coexistence.

## Principal uncertainties requiring calibration

1. **No direct human mucus corrinoid concentration** was found; fecal values cannot be assumed to equal mucus activity.
2. **No receptor-specific isogenic μmax measurements in colon-like medium** support exact BtuB-, FepA-, or double-null percentages.
3. **No immunity-gene-only competition experiment** supports a precise expression cost.
4. **No direct colicin diffusion measurement in native colonic mucus** was found.
5. **No reliable direct burst-size measurement in molecules per lysed cell** was found.
6. **No universal ColB conjugation coefficient** exists; plasmid architecture, helper elements, host, and spatial interface dominate.
7. **Healthy-human mucus-associated E. coli density and microcolony spacing** remain poorly quantified. Fecal density should not be presented as a mucosal measurement.

## Key bibliography (BibTeX)

```bibtex
@article{Feldgarden1999,
  author={Feldgarden, Michael and Riley, Margaret A.},
  title={The phenotypic and fitness effects of colicin resistance in Escherichia coli K-12},
  journal={Evolution}, year={1999}, volume={53}, pages={1019--1027},
  doi={10.1111/j.1558-5646.1999.tb04517.x}
}
@article{Feldgarden1998,
  author={Feldgarden, Michael and Riley, Margaret A.},
  title={High levels of colicin resistance in Escherichia coli},
  journal={Evolution}, year={1998}, volume={52}, pages={1270--1276},
  doi={10.1111/j.1558-5646.1998.tb02008.x}
}
@article{Kerr2002,
  author={Kerr, Benjamin and Riley, Margaret A. and Feldman, Marcus W. and Bohannan, Brendan J. M.},
  title={Local dispersal promotes biodiversity in a real-life game of rock-paper-scissors},
  journal={Nature}, year={2002}, volume={418}, pages={171--174},
  doi={10.1038/nature00823}
}
@article{Allen2008,
  author={Allen, Robert H. and Stabler, Sally P.},
  title={Identification and quantitation of cobalamin and cobalamin analogues in human feces},
  journal={American Journal of Clinical Nutrition}, year={2008}, volume={87}, pages={1324--1335},
  doi={10.1093/ajcn/87.5.1324}
}
@article{Kurisu2003,
  author={Kurisu, Genji and Zakharov, Stanislav D. and Zhalnina, Mariya V. and others},
  title={The structure of BtuB with bound colicin E3 R-domain implies a translocon},
  journal={Nature Structural Biology}, year={2003}, volume={10}, pages={948--954},
  doi={10.1038/nsb997}
}
@article{Cherezov2006,
  author={Cherezov, Vadim and Yamashita, Eiki and Liu, Weizhong and Zhalnina, Mariya V. and Cramer, William A.},
  title={In meso structure of the cobalamin transporter BtuB at 1.95 Å resolution},
  journal={Journal of Molecular Biology}, year={2006}, volume={364}, pages={716--734},
  doi={10.1016/j.jmb.2006.09.022}
}
@article{Housden2005,
  author={Housden, Nicholas G. and Loftus, Steven R. and Moore, Geoffrey R. and James, Richard and Kleanthous, Colin},
  title={Cell entry mechanism of enzymatic bacterial colicins: porin recruitment and thermodynamics of receptor binding},
  journal={Proceedings of the National Academy of Sciences USA}, year={2005}, volume={102}, pages={13849--13854},
  doi={10.1073/pnas.0503567102}
}
@article{Cavard1994,
  author={Cavard, Danièle},
  title={Rescue by vitamin B12 of Escherichia coli cells treated with colicins A and E allows measurement of the kinetics of colicin binding on BtuB},
  journal={FEMS Microbiology Letters}, year={1994}, volume={116}, pages={37--42},
  doi={10.1111/j.1574-6968.1994.tb06672.x}
}
@article{CohenKhait2021,
  author={Cohen-Khait, Ruth and Harmalkar, Ameya and Pham, Phuong and others},
  title={Colicin-mediated transport of DNA through the iron transporter FepA},
  journal={mBio}, year={2021}, volume={12}, number={5},
  doi={10.1128/mBio.01787-21}
}
@article{Sharp2019,
  author={Sharp, Connor and Boinett, Christine and Cain, Amy and others},
  title={O-antigen-dependent colicin insensitivity of uropathogenic Escherichia coli},
  journal={Journal of Bacteriology}, year={2019}, volume={201}, number={4},
  doi={10.1128/JB.00545-18}
}
@article{Martinson2020,
  author={Martinson, Jonathan N. V. and Walk, Seth T.},
  title={Escherichia coli residency in the gut of healthy human adults},
  journal={EcoSal Plus}, year={2020}, volume={9}, number={1},
  doi={10.1128/ecosalplus.ESP-0003-2020}
}
@article{Condamine2025,
  author={Condamine, Bénédicte and Morel-Journel, Thibaut and Tesson, Florian and others},
  title={Strain phylogroup and environmental constraints shape Escherichia coli dynamics and diversity over a 20-year human gut time series},
  journal={ISME Journal}, year={2025}, volume={19},
  doi={10.1093/ismejo/wrae245}
}
@article{FosterNyarko2022,
  author={Foster-Nyarko, Ebenezer and Pallen, Mark J.},
  title={The microbial ecology of Escherichia coli in the vertebrate gut},
  journal={FEMS Microbiology Reviews}, year={2022}, volume={46}, number={3},
  doi={10.1093/femsre/fuac008}
}
@article{Macfarlane2005,
  author={Macfarlane, Sandra and Woodmansey, Emma J. and Macfarlane, George T.},
  title={Colonization of mucin by human intestinal bacteria and establishment of biofilm communities in a two-stage continuous culture system},
  journal={Applied and Environmental Microbiology}, year={2005}, volume={71}, pages={7483--7492},
  doi={10.1128/AEM.71.11.7483-7492.2005}
}
@article{Poxton1997,
  author={Poxton, Ian R. and Brown, Robert and Sawyerr, A. and Ferguson, Anne},
  title={Mucosa-associated bacterial flora of the human colon},
  journal={Journal of Medical Microbiology}, year={1997}, volume={46}, pages={85--91},
  doi={10.1099/00222615-46-1-85}
}
@article{Saltzman1994,
  author={Saltzman, W. Mark and Radomsky, Michael L. and Whaley, Kevin J. and Cone, Richard A.},
  title={Antibody diffusion in human cervical mucus},
  journal={Biophysical Journal}, year={1994}, volume={66}, pages={508--515},
  doi={10.1016/S0006-3495(94)80802-1}
}
@article{Yildiz2015,
  author={Yildiz, Hasan M. and McKelvey, Craig A. and Marsac, Patrick J. and Carrier, Rebecca L.},
  title={Size selectivity of intestinal mucus to diffusing particulates is dependent on surface chemistry and exposure to lipids},
  journal={Journal of Drug Targeting}, year={2015}, volume={23}, pages={768--774},
  doi={10.3109/1061186X.2015.1086359}
}
@article{Macierzanka2019,
  author={Macierzanka, Adam and Mackie, Alan R. and Krupa, Lukasz},
  title={Permeability of the small intestinal mucus for physiologically relevant studies: impact of mucus location and ex vivo treatment},
  journal={Scientific Reports}, year={2019}, volume={9},
  doi={10.1038/s41598-019-53933-5}
}
@article{Spriewald2015,
  author={Spriewald, Stefanie and Glaser, Jana and Beutler, Markus and Koeppel, Martin B. and Stecher, Bärbel},
  title={Reporters for single-cell analysis of colicin Ib expression in Salmonella enterica serovar Typhimurium},
  journal={PLoS ONE}, year={2015}, volume={10}, pages={e0144647},
  doi={10.1371/journal.pone.0144647}
}
@article{Spriewald2020,
  author={Spriewald, Stefanie and Stadler, Eva and Hense, Burkhard A. and others},
  title={Evolutionary stabilization of cooperative toxin production through a bacterium-plasmid-phage interplay},
  journal={mBio}, year={2020}, volume={11}, number={4},
  doi={10.1128/mBio.00912-20}
}
@article{Jones2021,
  author={Jones, Emma C. and Uphoff, Stephan},
  title={Single-molecule imaging of LexA degradation in Escherichia coli elucidates regulatory mechanisms and heterogeneity of the SOS response},
  journal={Nature Microbiology}, year={2021}, volume={6}, pages={981--990},
  doi={10.1038/s41564-021-00930-y}
}
@article{BayramogluGuven2022,
  author={Bayramoglu-Güven, Bihter and Ghazaryan, Lusine and Toubiana, David and Gillor, Osnat},
  title={Colicin E2 expression in Escherichia coli biofilms: induction and regulation revisited},
  journal={Current Research in Microbial Sciences}, year={2022}, volume={3}, pages={100171},
  doi={10.1016/j.crmicr.2022.100171}
}
@article{Lerminiaux2025,
  author={Lerminiaux, Nicole A. and Kaufman, Jaycee M. and Schnell, Laura J. and others},
  title={Lysis of Escherichia coli by colicin Ib contributes to bacterial cross-feeding by releasing active beta-galactosidase},
  journal={ISME Journal}, year={2025}, volume={19},
  doi={10.1093/ismejo/wraf032}
}
@article{Zwaig1967,
  author={de Zwaig, Rosa Nagel and Luria, Salvador E.},
  title={Genetics and physiology of colicin-tolerant mutants of Escherichia coli},
  journal={Journal of Bacteriology}, year={1967}, volume={94}, pages={1112--1123},
  doi={10.1128/JB.94.4.1112-1123.1967}
}
@article{Shannon1967,
  author={Shannon, R. and Hedges, A. J.},
  title={Kinetics of lethal adsorption of colicin E2 by Escherichia coli},
  journal={Journal of Bacteriology}, year={1967}, volume={93}, pages={1353--1359},
  doi={10.1128/JB.93.4.1353-1359.1967}
}
@article{Reynolds1969,
  author={Reynolds, B. L. and Reeves, P. R.},
  title={Kinetics of adsorption of colicin CA42-E2 and reversal of its bactericidal activity},
  journal={Journal of Bacteriology}, year={1969}, volume={100}, pages={301--309},
  doi={10.1128/JB.100.1.301-309.1969}
}
@article{Smallwood2009,
  author={Smallwood, Chuck R. and Marco, Amparo Gala and Xiao, Qiaobin and others},
  title={Fluoresceination of FepA during colicin B killing: effects of temperature, toxin and TonB},
  journal={Molecular Microbiology}, year={2009}, volume={72}, pages={1171--1180},
  doi={10.1111/j.1365-2958.2009.06715.x}
}
@article{Cascales2007,
  author={Cascales, Eric and Buchanan, Susan K. and Duche, Denis and others},
  title={Colicin biology},
  journal={Microbiology and Molecular Biology Reviews}, year={2007}, volume={71}, pages={158--229},
  doi={10.1128/MMBR.00036-06}
}
@article{Riley1992,
  author={Riley, Margaret A. and Gordon, David M.},
  title={A survey of Col plasmids in natural isolates of Escherichia coli and an investigation into the stability of Col-plasmid lineages},
  journal={Journal of General Microbiology}, year={1992}, volume={138}, pages={1345--1352},
  doi={10.1099/00221287-138-7-1345}
}
@article{Johnson2006,
  author={Johnson, Timothy J. and Johnson, Sara J. and Nolan, Lisa K.},
  title={Complete DNA sequence of a ColBM plasmid from avian pathogenic Escherichia coli suggests that it evolved from closely related ColV virulence plasmids},
  journal={Journal of Bacteriology}, year={2006}, volume={188}, pages={5975--5983},
  doi={10.1128/JB.00204-06}
}
@article{PerezMendoza2009,
  author={Pérez-Mendoza, Daniel and de la Cruz, Fernando},
  title={Escherichia coli genes affecting recipient ability in plasmid conjugation: are there any?},
  journal={BMC Genomics}, year={2009}, volume={10}, pages={71},
  doi={10.1186/1471-2164-10-71}
}
@article{Element2023,
  author={Element, Sarah J. and Moran, Robert A. and Beattie, Emilie and others},
  title={Growth in a biofilm promotes conjugation of a blaNDM-1-bearing plasmid between Klebsiella pneumoniae strains},
  journal={mSphere}, year={2023}, volume={8}, number={4},
  doi={10.1128/msphere.00170-23}
}
@article{Djermoun2025,
  author={Djermoun, Sarah and Rode, Daniel K. H. and Jiménez-Siebert, Eva and others},
  title={Biofilm architecture determines the dissemination of conjugative plasmids},
  journal={Proceedings of the National Academy of Sciences USA}, year={2025}, volume={122},
  doi={10.1073/pnas.2417452122}
}
```

The most defensible baseline is therefore a model with **small, environment-dependent resistance and immunity costs; nanomolar true cobalamin embedded in approximately micromolar total corrinoids; near-water protein diffusion unless reversible binding is included; about 1% spontaneous toxin-expression heterogeneity; one productive translocation as the lethal hit; and no ColE1 transfer in the absence of a helper**. Burst size, mucosal density, mucus corrinoids, and conjugation should be designated calibration parameters rather than fixed literature constants.

References

1.: Stefanie Spriewald, Eva Stadler, Burkhard A. Hense, Philipp C. Münch, Alice C. McHardy, Anna S. Weiss, Nancy Obeng, Johannes Müller, and Bärbel Stecher. Evolutionary stabilization of cooperative toxin production through a bacterium-plasmid-phage interplay. Aug 2020. URL: https://doi.org/10.1128/mbio.00912-20, doi:10.1128/mbio.00912-20. This article has 17 citations and is from a domain leading peer-reviewed journal.

2.: Robert H Allen and Sally P Stabler. Identification and quantitation of cobalamin and cobalamin analogues in human feces. The American journal of clinical nutrition, 87 5:1324-35, May 2008. URL: https://doi.org/10.1093/ajcn/87.5.1324, doi:10.1093/ajcn/87.5.1324. This article has 233 citations.

3.: Robert H Allen and Sally P Stabler. Identification and quantitation of cobalamin and cobalamin analogues in human feces. The American journal of clinical nutrition, 87 5:1324-35, May 2008. URL: https://doi.org/10.1093/ajcn/87.5.1324, doi:10.1093/ajcn/87.5.1324. This article has 233 citations.

4.: Robert H Allen and Sally P Stabler. Identification and quantitation of cobalamin and cobalamin analogues in human feces. The American journal of clinical nutrition, 87 5:1324-35, May 2008. URL: https://doi.org/10.1093/ajcn/87.5.1324, doi:10.1093/ajcn/87.5.1324. This article has 233 citations.

5.: Bénédicte Condamine, Thibaut Morel-Journel, Florian Tesson, Guilhem Royer, Mélanie Magnan, Aude Bernheim, Erick Denamur, François Blanquart, and Olivier Clermont. Strain phylogroup and environmental constraints shape <i>escherichia coli</i> dynamics and diversity over a 20-year human gut time series. The ISME Journal, Dec 2025. URL: https://doi.org/10.1093/ismejo/wrae245, doi:10.1093/ismejo/wrae245. This article has 14 citations.

6.: Jonathan N. V. Martinson and Seth T. Walk. <i>escherichia coli</i> residency in the gut of healthy human adults. Dec 2020. URL: https://doi.org/10.1128/ecosalplus.esp-0003-2020, doi:10.1128/ecosalplus.esp-0003-2020. This article has 281 citations.

7.: Jonathan N. V. Martinson and Seth T. Walk. <i>escherichia coli</i> residency in the gut of healthy human adults. Dec 2020. URL: https://doi.org/10.1128/ecosalplus.esp-0003-2020, doi:10.1128/ecosalplus.esp-0003-2020. This article has 281 citations.

8.: Bénédicte Condamine, Thibaut Morel-Journel, Florian Tesson, Guilhem Royer, Mélanie Magnan, Aude Bernheim, Erick Denamur, François Blanquart, and Olivier Clermont. Strain phylogroup and environmental constraints shape <i>escherichia coli</i> dynamics and diversity over a 20-year human gut time series. The ISME Journal, Dec 2025. URL: https://doi.org/10.1093/ismejo/wrae245, doi:10.1093/ismejo/wrae245. This article has 14 citations.

9.: Jonathan N. V. Martinson and Seth T. Walk. <i>escherichia coli</i> residency in the gut of healthy human adults. Dec 2020. URL: https://doi.org/10.1128/ecosalplus.esp-0003-2020, doi:10.1128/ecosalplus.esp-0003-2020. This article has 281 citations.

10.: W.M. Saltzman, M.L. Radomsky, K.J. Whaley, and R.A. Cone. Antibody diffusion in human cervical mucus. Biophysical journal, 66 2 Pt 1:508-15, Feb 1994. URL: https://doi.org/10.1016/s0006-3495(94)80802-1, doi:10.1016/s0006-3495(94)80802-1. This article has 484 citations and is from a domain leading peer-reviewed journal.

11.: W.M. Saltzman, M.L. Radomsky, K.J. Whaley, and R.A. Cone. Antibody diffusion in human cervical mucus. Biophysical journal, 66 2 Pt 1:508-15, Feb 1994. URL: https://doi.org/10.1016/s0006-3495(94)80802-1, doi:10.1016/s0006-3495(94)80802-1. This article has 484 citations and is from a domain leading peer-reviewed journal.

12.: Stefanie Spriewald, Jana Glaser, Markus Beutler, Martin B. Koeppel, and Bärbel Stecher. Reporters for single-cell analysis of colicin ib expression in salmonella enterica serovar typhimurium. PLoS ONE, 10:e0144647, Dec 2015. URL: https://doi.org/10.1371/journal.pone.0144647, doi:10.1371/journal.pone.0144647. This article has 19 citations and is from a peer-reviewed journal.

13.: Nicole A Lerminiaux, Jaycee M Kaufman, Laura J Schnell, Sean D Workman, Danae M Suchan, Carsten Kröger, Brian P Ingalls, and Andrew D S Cameron. Lysis of <i>escherichia coli</i> by colicin ib contributes to bacterial cross-feeding by releasing active β-galactosidase. The ISME Journal, Feb 2025. URL: https://doi.org/10.1093/ismejo/wraf032, doi:10.1093/ismejo/wraf032. This article has 5 citations.

14.: W.M. Saltzman, M.L. Radomsky, K.J. Whaley, and R.A. Cone. Antibody diffusion in human cervical mucus. Biophysical journal, 66 2 Pt 1:508-15, Feb 1994. URL: https://doi.org/10.1016/s0006-3495(94)80802-1, doi:10.1016/s0006-3495(94)80802-1. This article has 484 citations and is from a domain leading peer-reviewed journal.