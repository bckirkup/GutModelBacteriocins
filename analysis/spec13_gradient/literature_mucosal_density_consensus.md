# I am parameterizing an individual-based model of Enterobacteriaceae in the human

**Question:** I am parameterizing an individual-based model of Enterobacteriaceae in the human colonic mucus layer. Two measured quantities needed, each with units and the method used to measure it: (1) The total bacterial density in the intestinal/colonic MUCUS layer or mucosa-associated compartment — CFU per mL or cells per mL of mucus (NOT luminal/fecal density). Relevant settings: healthy human colonic mucosa, mucus layer, biopsy-adjacent mucus, or mucus gel. Report the range across studies and which methods (qPCR, FISH counts, CFU plating, 16S sequencing normalized) produced which values. (2) The fraction of the total mucus/mucosa-associated bacterial population that is Enterobacteriaceae (or facultative anaerobes / Proteobacteria as reported), in healthy adults — as a proportion or percentage, with the measurement method. Also useful: whether Enterobacteriaceae fraction differs between mucus and stool/lumen in the same subjects. Report measured values only, with study names, n, and where in each paper the number appears. If the mucus-specific density is not directly measured in humans, say so and give the closest measured proxy (e.g., mucosa-associated counts per gram of tissue or per biopsy).

# Direct Measurements of **Colonic Mucus** Enterobacteriaceae Parameters

## 1. Introduction

For parameterizing an individual-based model of **Enterobacteriaceae in the human colonic mucus layer**, the literature supports one firm conclusion: healthy-human studies rarely report **native mucus bacterial density as cells or CFU per mL of mucus**, and most quantitative data instead use proxies such as **gene copies per mg of mucus**, **genomes per g of biopsy**, **cells per biopsy**, or **microscopic counts per mucus area or epithelial length** [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[2.1]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[3.1]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/) [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) [[5.1]](https://consensus.app/papers/comparative-study-of-the-intestinal-mucus-barrier-in-swidsinski-loening-baucke/be881b476d0e5409b29bdaa7d8f6faa9/). Healthy mucosa-associated communities are consistently distinct from stool or lumen, and across sequencing-based studies the mucosal compartment tends to contain a **higher Proteobacteria fraction** than stool, while Enterobacteriaceae-specific fractions are reported only in a few healthy cohorts [[6.1]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[7.1]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/) [[8.1]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[8.2]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[2.2]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/).

The most directly relevant human measurements come from colonic mucus gel brushings quantified by qPCR, classic quantitative culture of mucosal biopsies, and FISH-based microscopy of mucosal compartments [[1.2]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[4.2]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) [[5.2]](https://consensus.app/papers/comparative-study-of-the-intestinal-mucus-barrier-in-swidsinski-loening-baucke/be881b476d0e5409b29bdaa7d8f6faa9/). These methods are not interchangeable: qPCR yields 16S copy number or inferred genomes, culture yields viable CFU, FISH yields in situ cell counts or fractions, and 16S sequencing yields relative abundance rather than density [[2.3]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[4.3]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) [[5.2]](https://consensus.app/papers/comparative-study-of-the-intestinal-mucus-barrier-in-swidsinski-loening-baucke/be881b476d0e5409b29bdaa7d8f6faa9/) [[9.1]](https://consensus.app/papers/composition-of-the-mucosaassociated-microbiota-along-the-vuik-dicksved/5c08f2c41cae58a5bd448cf58d17bad5/). For modeling, the safest approach is to separate **direct mucus measurements** from **mucosal-biopsy proxies**, and to keep **Proteobacteria-level fractions** distinct from **Enterobacteriaceae-specific measurements** [[10.1]](https://consensus.app/papers/depthdependent-differences-in-community-structure-of-the-lavelle-lennon/a1a7d0067f51571e84c9fbafe4666d58/) [[11.1]](https://consensus.app/papers/analysis-of-colonic-mucosaassociated-microbiota-using-miyauchi-taida/68883bb0605a5302ac56f914395ffb27/) [[12.1]](https://consensus.app/papers/relationship-between-mucosaassociated-gut-microbiota-juge/22f5cb8ea81151b9ba5acdaa649eb25a/).


**Figure 1:** Consensus on direct human mucus parameter availability

Direct measurements exist, but they are sparse and mostly indirect for absolute density. Enterobacteriaceae fractions are even less directly quantified than total mucosal bacterial burden.

## 2. Methods

The search ran over 220 million research papers in Consensus, spanning Semantic Scholar, PubMed, and other sources. The corpus was screened for healthy human colonic mucus, mucosa-associated microbiota, biopsy-adjacent mucus, brushings, lavages, and paired stool–mucosa comparisons, with emphasis on studies reporting quantitative units and explicit measurement methods. Extracted values were limited to measured numbers in healthy adults or healthy control groups, and were separated into native mucus measurements versus tissue or biopsy proxies.

| Retrieved | Eligible | Included |
|-----------|----------|----------|
| 297.9K | 212 | 100 | 

**Figure 2:** Search funnel for mucus-associated microbiota studies

Searches combined mucus-layer density terms, Enterobacteriaceae/Proteobacteria abundance terms, paired stool–mucosa comparisons, and citation-following from quantitative benchmark studies.

## 3. Results

### 3.1 Key Papers

A small set of papers anchor the quantitative literature because they either provide the closest direct human mucus measurements or the clearest healthy mucosa-versus-stool comparisons [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[6.1]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/). No single paper provides both requested quantities in the same native-mucus unit system, so triangulation across methods is required [[12.1]](https://consensus.app/papers/relationship-between-mucosaassociated-gut-microbiota-juge/22f5cb8ea81151b9ba5acdaa649eb25a/) [[12.2]](https://consensus.app/papers/relationship-between-mucosaassociated-gut-microbiota-juge/22f5cb8ea81151b9ba5acdaa649eb25a/).

| Paper | Summary |
|---|---|
| [[1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) | Mucus qPCR copies per mg [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) |
| [[6]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) | Mucosa vs stool Proteobacteria fractions [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) |
| [[4]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) | Viable CFU per biopsy proxy [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) |

**Figure 3:** Foundational quantitative papers for requested parameters

### 3.2 Absolute Density

Healthy-human **native colonic mucus** was directly quantified by RT-qPCR in mucus gel brushings in Earley et al., with total bacterial 16S copy number reported per mg of mucus rather than per mL; for healthy controls, median **R. hominis** copy number was 1.13×10^5 per mg after normalization context to total bacterial copy number, and the same study states species abundance was normalized to total bacterial copy number in mucus gel brushings from 20 healthy controls, 20 active UC, and 14 quiescent UC subjects [[1.2]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[1.3]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/). This is the closest direct healthy-human mucus-layer denominator in the supplied corpus, but it is **copies per mg mucus**, not cells per mL mucus [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/).

Older healthy-human mucosal culture studies instead reported **per-biopsy** counts. Poxton et al. studied 6 non-inflammatory controls and found viable aerobic counts of **2.4×10^3 to 1.3×10^6 CFU per biopsy sample** and total anaerobic counts of **1.4×10^5 to 3×10^7 CFU per biopsy sample**, with each biopsy weighing about 5.6 mg [[4.2]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/). Zoetendal et al. reported **10^5 to 10^6 bacteria per biopsy sample** by flow cytometry across colon biopsy specimens from 10 individuals, again as a biopsy-level proxy rather than mucus volume density [[3.2]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/) [[3.1]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/).

Lyra et al. provide another tissue proxy using 16S-targeted qPCR expressed as **log10 genomes per g wet sample**. In endoscopically and histologically normal mucosal biopsies, total eubacteria were quantified by qPCR with results expressed as log10 genomes per gram after host-DNA correction [[2.1]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[2.4]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[2.5]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/). Swidsinski’s FISH work counted bacteria in a defined microscopic mucus square or reported mean mucosal-bacteria concentration in units of **10^10/mL**, but these refer to microscopic mucosal biofilm or adherent mucus compartments rather than recovered native mucus volume [[5.1]](https://consensus.app/papers/comparative-study-of-the-intestinal-mucus-barrier-in-swidsinski-loening-baucke/be881b476d0e5409b29bdaa7d8f6faa9/) [[13.1]](https://consensus.app/papers/spatial-organization-and-composition-of-the-mucosal-flora-swidsinski-weber/d198671ffa7052e985671f7ba5b6d38e/) [[13.2]](https://consensus.app/papers/spatial-organization-and-composition-of-the-mucosal-flora-swidsinski-weber/d198671ffa7052e985671f7ba5b6d38e/).

**Practical range for quantity (1), keeping denominators separate**

| Study | Healthy setting, n | Method | Reported value | Where in paper |
|---|---|---|---|---|
| Earley 2020 | Colonic mucus gel brushings, HC n=20 | RT-qPCR | total bacterial copy number per mg mucus; mucus-denominator study, not per mL [[1.2]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) | Table 1 / results table context |
| Poxton 1997 | Proximal colon and rectal biopsies, non-inflammatory n=6 | Quantitative culture | aerobes 2.4×10^3–1.3×10^6 CFU/biopsy; anaerobes 1.4×10^5–3×10^7 CFU/biopsy [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) | Abstract results |
| Zoetendal 2002 | Colon biopsies, n=10 | Flow cytometry | 10^5–10^6 bacteria per biopsy sample [[3.1]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/) | Abstract results |
| Lyra 2012 | Right and left colon biopsies, n=20 | qPCR | eubacteria as log10 genomes/g wet sample [[2.1]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[2.6]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) | Methods, Fig. 1/Table 3 context |
| Swidsinski 2005/2006 | Healthy controls | FISH microscopy | mucus/adherent bacteria in microscopic compartments; mean control values around 0.02×10^10/mL in controls in one dataset [[13.1]](https://consensus.app/papers/spatial-organization-and-composition-of-the-mucosal-flora-swidsinski-weber/d198671ffa7052e985671f7ba5b6d38e/) [[13.2]](https://consensus.app/papers/spatial-organization-and-composition-of-the-mucosal-flora-swidsinski-weber/d198671ffa7052e985671f7ba5b6d38e/) | Table 4 / abstract |

**Figure 4:** Absolute density estimates use incompatible denominators

### 3.3 Enterobacteriaceae and Proteobacteria Fraction

Healthy-adult **Enterobacteriaceae-specific** percentages are reported in Hong et al. using 16S pyrosequencing of 30 biopsies from 5 healthy individuals. Their Table 2 gives **unclassified Enterobacteriaceae** mean relative abundances by individual of **3.52%, 0.68%, 3.37%, 1.39%, and 3.81%**, implying a cross-individual span of **0.68% to 3.81%** in healthy mucosal biopsies [[8.3]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[8.2]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/). The same cohort had overall **Proteobacteria 8.6% ± 4.7%** of the mucosal community, with individual biopsies ranging **1.7% to 19.2%** [[8.1]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[8.4]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/).

Several healthy cohorts report **Proteobacteria**, rather than Enterobacteriaceae, as the usable facultative-anaerobe proxy. Ringel et al. found Proteobacteria at **19.3% in distal colonic mucosa-associated microbiota versus 9.3% in feces** in 24 healthy subjects using 16S sequencing of unprepped distal colonic biopsies and same-day feces [[6.3]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/). Borgo et al. likewise found sigmoid brush-sample MAM enriched in aerotolerant Proteobacteria relative to feces in 40 adults [[7.2]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/) [[7.1]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/). Hong et al. and Shen et al. reported mucosal Proteobacteria means of **8.6%** and **11%** respectively in broadly healthy or control mucosal biopsies [[8.1]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[14.1]](https://consensus.app/papers/molecular-characterization-of-mucosal-adherent-bacteria-shen-rawls/87f2030792715a7a83597700022a3129/).

Targeted qPCR gives an Enterobacteriaceae proxy on mucosal biopsies in healthy adults. Lyra et al. quantified Enterobacteriaceae at **4.68 ± 1.21 log10 genomes/g** in right-colon mucosa, **4.71 ± 0.83 log10 genomes/g** in left-colon mucosa, and **5.70 ± 2.00 log10 genomes/g** in feces; fecal levels trended higher, but the difference was not statistically significant for Enterobacteriaceae alone [[2.2]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[2.7]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[2.8]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/). E. coli itself was detected in fewer than half of mucosal samples but in all fecal samples in that study [[2.9]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[2.10]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/).

| Dimension | Healthy mucus/mucosa value | Stool/lumen comparator | Method | Citations |
|---|---|---|---|---|
| Enterobacteriaceae fraction | 0.68%–3.81% across 5 individuals | no paired stool in same subjects | 16S pyrosequencing biopsy relative abundance | [[8.2]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) |
| Proteobacteria fraction | 8.6% ± 4.7%; 1.7%–19.2% across biopsies | pooled stool outgroup only | 16S pyrosequencing biopsy relative abundance | [[8.1]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[8.4]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) |
| Proteobacteria fraction | 19.3% mucosa | 9.3% feces | 16S sequencing | [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) |
| Proteobacteria enrichment | MAM enriched | feces depleted | 16S sequencing brush vs stool | [[7.2]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/) [[7.1]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/) |
| Enterobacteriaceae quantity | 4.68–4.71 log10 genomes/g mucosa | 5.70 log10 genomes/g feces | targeted qPCR | [[2.2]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) |

**Figure 5:** Healthy mucosa shows modest Enterobacteriaceae and higher Proteobacteria

### 3.4 Additional Comparative Signals

Multiple studies support the broader point that mucosa-associated communities differ from stool and often show relative enrichment of facultative anaerobes or Proteobacteria in the mucosal compartment [[11.2]](https://consensus.app/papers/analysis-of-colonic-mucosaassociated-microbiota-using-miyauchi-taida/68883bb0605a5302ac56f914395ffb27/) [[8.5]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[12.2]](https://consensus.app/papers/relationship-between-mucosaassociated-gut-microbiota-juge/22f5cb8ea81151b9ba5acdaa649eb25a/) [[15.1]](https://consensus.app/papers/comparison-of-microbial-signatures-between-paired-faecal-mukhopadhya-martin/b0270458194856efa3d5ebafb8d17719/). Yilmaz et al. add spatial resolution within the mucosal interface, showing that aspirated superficial mucus in healthy volunteers was enriched for Enterobacteriaceae relative to deeper biopsy-associated communities, indicating that any single “mucus Enterobacteriaceae fraction” depends on whether one samples outer mucus, deeper mucus, or biopsy-adherent mucosa [[16.1]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/) [[16.2]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/).

Other healthy-control sequencing cohorts reinforce that Proteobacteria are a minority but recurrent mucosal phylum. Hong et al. reported 8.6% mucosal Proteobacteria, Shen et al. 11%, and Ringel et al. 19.3% in distal mucosa, while Flynn et al. found proximal mucosa relatively enriched in Enterobacteriaceae compared with distal mucosa in unprepped healthy colons [[8.1]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[14.1]](https://consensus.app/papers/molecular-characterization-of-mucosal-adherent-bacteria-shen-rawls/87f2030792715a7a83597700022a3129/) [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[17.1]](https://consensus.app/papers/spatial-variation-of-the-native-colon-microbiota-in-flynn-ruffin/bbbbc69051a65db1865fd6efb4b1e230/). de Cárcer et al. similarly reported that several Enterobacteriaceae increase in relative abundance toward the rectum [[18.1]](https://consensus.app/papers/numerical-ecology-validates-a-biogeographical-cárcer-cuív/4490e0e8ea9a539f961bb76942a21f82/).

### Results Timeline

Colonic mucus-associated quantification has moved from biopsy culture and microscopy to spatially resolved sequencing and targeted qPCR.

- **1997**
  - 1 paper: [[4]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/)
- **2002**
  - 1 paper: [[3]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/)
- **2005**
  - 2 papers: [[13]](https://consensus.app/papers/spatial-organization-and-composition-of-the-mucosal-flora-swidsinski-weber/d198671ffa7052e985671f7ba5b6d38e/) [[19]](https://consensus.app/papers/molecular-characterization-of-rectal-mucosa‐associated-mylonaki-rayment/0c991d06ce8651eba45fb7e892deeb93/)
- **2006**
  - 1 paper: [[5]](https://consensus.app/papers/comparative-study-of-the-intestinal-mucus-barrier-in-swidsinski-loening-baucke/be881b476d0e5409b29bdaa7d8f6faa9/)
- **2010**
  - 1 paper: [[14]](https://consensus.app/papers/molecular-characterization-of-mucosal-adherent-bacteria-shen-rawls/87f2030792715a7a83597700022a3129/)
- **2011**
  - 2 papers: [[8]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[18]](https://consensus.app/papers/numerical-ecology-validates-a-biogeographical-cárcer-cuív/4490e0e8ea9a539f961bb76942a21f82/)
- **2012**
  - 1 paper: [[2]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/)
- **2013**
  - 1 paper: [[10]](https://consensus.app/papers/depthdependent-differences-in-community-structure-of-the-lavelle-lennon/a1a7d0067f51571e84c9fbafe4666d58/)
- **2015**
  - 1 paper: [[6]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/)
- **2017**
  - 1 paper: [[17]](https://consensus.app/papers/spatial-variation-of-the-native-colon-microbiota-in-flynn-ruffin/bbbbc69051a65db1865fd6efb4b1e230/)
- **2018**
  - 1 paper: [[7]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/)
- **2019**
  - 2 papers: [[9]](https://consensus.app/papers/composition-of-the-mucosaassociated-microbiota-along-the-vuik-dicksved/5c08f2c41cae58a5bd448cf58d17bad5/) [[20]](https://consensus.app/papers/analysis-of-the-colonic-mucosa-associated-microbiota-mam-matsumoto-kuroki/a06321fda7bf5e5badc218f217793d03/)
- **2020**
  - 1 paper: [[1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/)
- **2022**
  - 3 papers: [[11]](https://consensus.app/papers/analysis-of-colonic-mucosaassociated-microbiota-using-miyauchi-taida/68883bb0605a5302ac56f914395ffb27/) [[12]](https://consensus.app/papers/relationship-between-mucosaassociated-gut-microbiota-juge/22f5cb8ea81151b9ba5acdaa649eb25a/) [[15]](https://consensus.app/papers/comparison-of-microbial-signatures-between-paired-faecal-mukhopadhya-martin/b0270458194856efa3d5ebafb8d17719/)
- **2026**
  - 1 paper: [[16]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/)
**Figure 6:** Timeline of mucus-associated microbiota quantification studies

The arc is methodological rather than purely biological: early work quantified viable or microscopic burdens, mid-period work established mucosa–stool divergence, and newer work resolved superficial mucus versus deeper mucosal niches. Larger markers indicate more citations.

### Top Contributors

| Type | Name | Papers |
|------|------|--------|
| Author | A. Lavelle | [[21]](https://consensus.app/papers/spatial-variation-of-the-colonic-microbiota-in-patients-lavelle-lennon/cd3a62c8bd7e5aaba9162929b9edfe4d/) [[10]](https://consensus.app/papers/depthdependent-differences-in-community-structure-of-the-lavelle-lennon/a1a7d0067f51571e84c9fbafe4666d58/) [[22]](https://consensus.app/papers/correlations-between-colonic-crypt-mucin-chemotype-lennon-bambury/c795a176467a5aebbe423b9165f14a40/) |
| Author | A. Swidsinski | [[5]](https://consensus.app/papers/comparative-study-of-the-intestinal-mucus-barrier-in-swidsinski-loening-baucke/be881b476d0e5409b29bdaa7d8f6faa9/) [[13]](https://consensus.app/papers/spatial-organization-and-composition-of-the-mucosal-flora-swidsinski-weber/d198671ffa7052e985671f7ba5b6d38e/) [[23]](https://consensus.app/papers/mucosal-flora-in-inflammatory-bowel-disease-swidsinski-ladhoff/475fe520b86253179e7eeb2473c17a02/) |
| Author | S. Macfarlane | [[24]](https://consensus.app/papers/mucosaassociated-bacterial-diversity-in-relation-to-ahmed-macfarlane/0a2f25bc1d5357c58234fefcd36c767c/) [[25]](https://consensus.app/papers/colonization-of-mucin-by-human-intestinal-bacteria-and-macfarlane-woodmansey/5ed5d514addf50b8adf2ade1d25e4ac6/) |
| Journal | *Gut* | [[21]](https://consensus.app/papers/spatial-variation-of-the-colonic-microbiota-in-patients-lavelle-lennon/cd3a62c8bd7e5aaba9162929b9edfe4d/) [[5]](https://consensus.app/papers/comparative-study-of-the-intestinal-mucus-barrier-in-swidsinski-loening-baucke/be881b476d0e5409b29bdaa7d8f6faa9/) [[26]](https://consensus.app/papers/thickness-of-adherent-mucus-gel-on-colonic-mucosa-in-humans-pullan-thomas/b04f01c1dc3c539d9d24e32d94a4481c/) |
| Journal | *Applied and Environmental Microbiology* | [[24]](https://consensus.app/papers/mucosaassociated-bacterial-diversity-in-relation-to-ahmed-macfarlane/0a2f25bc1d5357c58234fefcd36c767c/) [[3]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/) [[25]](https://consensus.app/papers/colonization-of-mucin-by-human-intestinal-bacteria-and-macfarlane-woodmansey/5ed5d514addf50b8adf2ade1d25e4ac6/) |
| Journal | *Gut Microbes* | [[16]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/) [[6]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[27]](https://consensus.app/papers/mucosal-washes-are-useful-for-sampling-intestinal-martinez-medina-ghazisaeedi/3bdf0537a60655f28727ebb3e1a7c1f7/) |

**Figure 7:** Authors and journals appearing most often

## 4. Discussion

The evidence base is strongest for the qualitative conclusion that healthy colonic mucosa or mucus harbors a microbiota distinct from stool and relatively enriched in Proteobacteria or other aerotolerant taxa [[11.2]](https://consensus.app/papers/analysis-of-colonic-mucosaassociated-microbiota-using-miyauchi-taida/68883bb0605a5302ac56f914395ffb27/) [[6.4]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[7.1]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/) [[12.1]](https://consensus.app/papers/relationship-between-mucosaassociated-gut-microbiota-juge/22f5cb8ea81151b9ba5acdaa649eb25a/). It is moderate for quantitative estimates of **Enterobacteriaceae or Proteobacteria fraction**, because several healthy cohorts report convergent relative-abundance values, but they use different sampling compartments and sequencing platforms [[8.2]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[14.1]](https://consensus.app/papers/molecular-characterization-of-mucosal-adherent-bacteria-shen-rawls/87f2030792715a7a83597700022a3129/) [[17.1]](https://consensus.app/papers/spatial-variation-of-the-native-colon-microbiota-in-flynn-ruffin/bbbbc69051a65db1865fd6efb4b1e230/).

The evidence is weaker for **absolute bacterial density in native human colonic mucus volume units**. The literature does provide measured burdens, but denominators differ sharply: CFU per biopsy, bacteria per biopsy, genomes per gram tissue, copy number per mg mucus, and microscopic counts per unit epithelial length or imaged mucus area [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) [[3.1]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/) [[2.1]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[19.1]](https://consensus.app/papers/molecular-characterization-of-rectal-mucosa‐associated-mylonaki-rayment/0c991d06ce8651eba45fb7e892deeb93/). That makes direct conversion to cells per mL mucus underdetermined without an external assumption about mucus density, recovery efficiency, biopsy geometry, 16S copy number per genome, and viable fraction.

A further limitation is compartment mixing. Biopsies include epithelium and some submucosa, brushes sample superficial mucus, lavages are intermediate, and recent spatial work indicates Enterobacteriaceae are concentrated mainly in superficial mucus rather than deeper mucosal surface communities [[10.1]](https://consensus.app/papers/depthdependent-differences-in-community-structure-of-the-lavelle-lennon/a1a7d0067f51571e84c9fbafe4666d58/) [[11.3]](https://consensus.app/papers/analysis-of-colonic-mucosaassociated-microbiota-using-miyauchi-taida/68883bb0605a5302ac56f914395ffb27/) [[16.1]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/) [[16.2]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/). For modeling, this means the “right” parameter depends on whether the simulated niche is **outer mucus**, **whole mucus gel**, or **biopsy-associated mucosa**.

The practical takeaway is that the best directly supported healthy-human inputs are a mucus-gel qPCR burden in **copies per mg mucus** and a mucosal Enterobacteriaceae/Proteobacteria fraction in the low-single-digit to low-double-digit percent range, with stool usually lower for Proteobacteria than mucosa [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[8.2]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/).

| Claim | Evidence Strength | Reasoning | Papers |
|---|---|---|---|
| Healthy human studies rarely report native colonic mucus density as cells or CFU per mL | Evidence strength: Strong (8/10) | Multiple papers quantify mucus-associated burden, but almost always with non-volume denominators | [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[2.1]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[3.1]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/) [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) |
| Closest direct absolute mucus measurement is qPCR copy number per mg mucus | Evidence strength: Moderate (7/10) | Direct mucus brushings exist, but reported in copies/mg rather than cells/mL | [[1.2]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) |
| Healthy mucosa has modest Enterobacteriaceae and higher Proteobacteria than stool | Evidence strength: Moderate (7/10) | Several healthy cohorts agree on Proteobacteria enrichment; fewer give Enterobacteriaceae specifically | [[8.2]](https://consensus.app/papers/pyrosequencingbased-analysis-of-the-mucosal-microbiota-hong-croix/3109553bb9ab5e9cad45ffaf204cdd17/) [[6.2]](https://consensus.app/papers/high-throughput-sequencing-reveals-distinct-microbial-ringel-maharshak/33e0af0f1efc55d2b0cf447374e75c0f/) [[7.1]](https://consensus.app/papers/body-mass-index-and-sex-affect-diverse-microbial-niches-borgo-garbossa/011f4b3732ac54b4858c610b6b278c5a/) |
| Superficial mucus is richer in Enterobacteriaceae than deeper mucosa | Evidence strength: Moderate (6/10) | Spatially resolved human evidence is recent and small but directly on-point | [[16.1]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/) [[16.2]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/) |
| A single universal conversion from biopsy counts to cells/mL mucus is currently justified | Evidence strength: Weak (2/10) | Denominator incompatibility and compartment mixing prevent a defensible direct conversion | [[10.1]](https://consensus.app/papers/depthdependent-differences-in-community-structure-of-the-lavelle-lennon/a1a7d0067f51571e84c9fbafe4666d58/) [[11.1]](https://consensus.app/papers/analysis-of-colonic-mucosaassociated-microbiota-using-miyauchi-taida/68883bb0605a5302ac56f914395ffb27/) [[12.2]](https://consensus.app/papers/relationship-between-mucosaassociated-gut-microbiota-juge/22f5cb8ea81151b9ba5acdaa649eb25a/) |

**Figure 8:** Key claims and evidence across included papers

## 5. Conclusion

For the two requested quantities, the literature gives a usable answer for **fractional abundance** and only a proxy-level answer for **absolute density**. Healthy-human native colonic mucus has been quantified directly, but mainly as **qPCR copy number per mg mucus**, while most other absolute measurements are per biopsy, per gram tissue, or microscopic compartment counts [[1.1]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) [[3.1]](https://consensus.app/papers/mucosaassociated-bacteria-in-the-human-gastrointestinal-zoetendal-wright/8d5ce928dc155223b3cbd5993a8e95f3/).

A cautious parameter set would therefore distinguish: **(i)** direct mucus burden in copies/mg mucus, **(ii)** biopsy-level or tissue-mass absolute proxies, and **(iii)** Enterobacteriaceae/Proteobacteria fraction by sampling compartment.

### Research Gaps

The coverage landscape is dominated by compositional studies and only sparsely populated by absolute quantification studies. The most important missing cell is a healthy-human study that measures **total cells per mL native colonic mucus and Enterobacteriaceae cells per mL in the same samples**.

| | Absolute Load | Relative Fraction | Spatial Depth | Healthy Adults | Enterobacteria |
|---|---|---|---|---|---|
| Brushings |  |  |  |  |  |
| Biopsies |  |  |  |  |  |
| Lavages |  |  |  |  |  |
| Mucus Gel |  |  |  |  |  |
| Paired Stool |  |  |  |  |  |

The sharpest gap is the absence of a common denominator across methods. Brush and mucus-gel studies resolve the relevant niche but usually report gene copies or relative abundance rather than cells per mL [[1.2]](https://consensus.app/papers/colonisation-of-the-colonic-mucus-gel-layer-with-earley-lennon/54d2eed44f85504c84609d968333a44b/) [[20.1]](https://consensus.app/papers/analysis-of-the-colonic-mucosa-associated-microbiota-mam-matsumoto-kuroki/a06321fda7bf5e5badc218f217793d03/) [[11.4]](https://consensus.app/papers/analysis-of-colonic-mucosaassociated-microbiota-using-miyauchi-taida/68883bb0605a5302ac56f914395ffb27/). Biopsy studies provide absolute burden but mix mucus with host tissue and therefore do not map cleanly onto a mucus-volume model state variable [[4.1]](https://consensus.app/papers/mucosaassociated-bacterial-flora-of-the-human-colon-poxton-brown/187a9195ef635752802036494f9cdbd0/) [[2.1]](https://consensus.app/papers/comparison-of-bacterial-quantities-in-left-and-right-colon-lyra-forssten/5fed58daf9635bbf8214066f2b0ec5e7/) [[10.1]](https://consensus.app/papers/depthdependent-differences-in-community-structure-of-the-lavelle-lennon/a1a7d0067f51571e84c9fbafe4666d58/). Spatial-depth resolution is also recent and limited, even though Yilmaz et al. suggest that superficial mucus and deeper mucosa differ materially in Enterobacteriaceae burden [[16.1]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/) [[16.2]](https://consensus.app/papers/spatial-mapping-of-human-colonic-niches-reveals-rapid-yilmaz-moulin/ffec438d399a5ac9b79935bd4d19b1c2/).

### Open Research Questions

The next advances are straightforward technically but have not yet been combined in the same healthy-human design.

| Question | Why |
|---|---|
| **What is the total bacterial cell density per mL of native healthy human colonic outer mucus?** | This would replace current copy-per-mg and per-biopsy proxies with the model’s needed state variable. |
| **What fraction of those mucus-resident cells are Enterobacteriaceae across outer versus deeper mucus layers?** | Spatially resolved fractions are needed because superficial mucus appears more Enterobacteriaceae-rich than deeper mucosa. |
| **How do qPCR, FISH, and viable culture calibrate against each other on the same human mucus sample?** | Cross-method calibration would allow robust conversion between total burden, viable burden, and taxon-specific fractions. |

For modeling **Enterobacteriaceae in the human colonic mucus layer**, the best-supported measured inputs are currently **copies per mg mucus** for total burden and **~0.7–3.8% Enterobacteriaceae or ~8.6–19.3% Proteobacteria** for healthy mucosa-associated relative abundance, with mucosa generally more Proteobacteria-rich than stool.
 
_These search results were found and analyzed using Consensus, an AI-powered search engine for research. Try it at https://consensus.app. © 2026 Consensus NLP, Inc. Personal, non-commercial use only; redistribution requires copyright holders’ consent._
 
## References
 
`1.` Earley H, Lennon G, Coffey JC, Winter D, O'Connell P. Colonisation of the colonic mucus gel layer with butyrogenic and hydrogenotropic bacteria in health and ulcerative colitis. *Scientific Reports.* 2020;11. doi:10.1038/s41598-021-86166-6
 
> [1.1] Table 3 Median copy number/mg and median relative abundance of R. hominis after normalisation for total bacterial copy number in each patient cohort. (Results)
 
> [1.2] Mucosal brushings were obtained from 20 healthy controls (HC), 20 patients with active colitis (AC) and 14 with quiescent colitis (QUC). (Abstract)
 
> [1.3] Abundance of each species was determined by RT-PCR. (Abstract)
 
`2.` Lyra A, Forssten S, Rolny P, et al. Comparison of bacterial quantities in left and right colon biopsies and faeces.. *World journal of gastroenterology.* 2012;18 32:4404-11. doi:10.3748/wjg.v18.i32.4404
 
> [2.1] Results were expressed as log 10 genomes per gram of sample (wet weight), taking into account the size and the 16S rDNA copy number of the standard species genome. (Methods)
 
> [2.2] ; 4.68 ± 1.21 log10/g, 4.71 ± 0.83 log10/g, 5.70 ± 2.00 log 10/g, P1 = 0.1927, P2 = 0.0605, P3 = 0.6476 for Enterobacteriaceae). (Results)
 
> [2.3] The qPCR reactions were performed using Ap - plied Biosystems Real-Time PCR system equipment (7500 Fast, Applied Biosystems, Foster City, CA, United States) and software applying in-house optimized assay conditions for the primer sequences presented in Table 2. (Methods)
 
> [2.4] Statistical analysis For mucosal samples, the proportion of host DNA was estimated according to the GAPDH qPCR result and subtracted prior to calculations. (Methods)
 
> [2.5] , Faecalibacte- rium prausnitzii (F. prausnitzii) and putative pathogens such as Escherichia coli (E. coli), Clostridium difficile (C. difficile), Helicobacter pylori (H. pylori ) and Staphylo- coccus aureus (S. aureus ) were analysed by quantita - tive polymerase chain reaction (qPCR). (Results)
 
> [2.6] Innovations and breakthroughs Real-time quantitative polymerase chain reaction (qPCR) allows independent Table 3 Pearson correlations between sample types Bacterial group/study period Right colon vs left colon Right colon vs faecal sample Left colon vs faecal sample Correlation coefficient P value Correlation coefficient P value Correlation coefficient P value Akkermansia muciniphila 0.14 0.63 -0.01 0.97 0.36 0.26 Bacteroidetes 0.61 0.02 0.45 0.12 0.17 0.61 Bifidobacterium 0.71 0.01 0.62 0.04 0.81 0.00 Clostridium Cluster Ⅳ 0.26 0.39 0.26 0.44 0.17 0.64 Clostridium Cluster ⅩⅣab 0.71 0.00 0.54 0.06 0.50 0.09 Collinsella aerofaciens 0.38 0.25 0.63 0.13 -0.87 0.13 Eubacteria 0.19 0.52 0.01 0.97 -0.31 0.33 Enterobacteriaceae 0.38 0.20 0.59 0.03 0.31 0.35 Faecalibacterium prausnitzii 0.55 0.04 0.76 0.00 0.28 0.38 Veillonella 0.33 0.46 0.64 0.09 0.46 0.54 Significant correlations (P < 0.05) are designated with bold font. (Results)
 
> [2.7] The bac - terial quantities between the two mucosal samples did not differ according to paired t-tests, whereas the faecal quantities of all analysed bacteria were sig - nificantly higher than those detected for either mucosal site ( P < 0.05), except for Enterobacteriaceae. (Results)
 
> [2.8] For Enterobacteriaceae the higher abundance in faeces was not statistically signifi - cant, but a similar trend was evident between the left side mucosal and faecal samples. (Discussion)
 
> [2.9] E. coli was detected in less than half of the mucosal samples at both sites, but was present in all faecal samples. (Results)
 
> [2.10] Escherichia coli (E. coli) was detected in less than half of the mucosal samples at both sites, while present in all faecal samples at log10 5.92 ± 1.04 genomes per gram of faeces. (Results)
 
`3.` Zoetendal E, Von Wright A, Vilpponen-Salmela T, Ben-Amor K, Akkermans A, De Vos WD. Mucosa-Associated Bacteria in the Human Gastrointestinal Tract Are Uniformly Distributed along the Colon and Differ from the Community Recovered from Feces. *Applied and Environmental Microbiology.* 2002;68:3401 - 3407. doi:10.1128/aem.68.7.3401-3407.2002
 
> [3.1] Flow cytometric analysis indicated that 10(5) to 10(6) bacteria were present in the biopsy samples. (Abstract)
 
> [3.2] TABLE 1. (Results)
 
`4.` Poxton I, Brown RH, Sawyerr AM, Ferguson A. Mucosa-associated bacterial flora of the human colon.. *Journal of medical microbiology.* 1997;46 1:85-91. doi:10.1099/00222615-46-1-85
 
> [4.1] Viable aerobic counts (aerobes and facultative organisms) ranged from 2.4 x 10(3) to 1.3 x 10(6) cfu/sample biopsy (5.6 mg) and total anaerobic counts were 10-10(2) times higher at (1.4 x 10(5))-(3 x 10(7)) cfu/sample. (Abstract)
 
> [4.2] Biopsy samples of mucosa were taken during colonoscopy from the proximal colon and rectum of 12 patients, six with ulcerative colitis (UC) and six with non-inflammatory conditions. (Abstract)
 
> [4.3] After anaerobic transport to the laboratory, biopsy specimens were examined by quantitative bacteriological culture on selective and non-selective media for total aerobic count, total anaerobic count, Bacteroides spp. (Abstract)
 
`5.` Swidsinski A, Loening-Baucke V, Theissig F, et al. Comparative study of the intestinal mucus barrier in normal and inflamed colon. *Gut.* 2006;56:343 - 350. doi:10.1136/gut.2006.098160
 
> [5.1] 8. Concentrations of bacteria within mucus: numbers of mucus-penetrating bacteria were calculated within a square of 10 6 10 m m, which was placed over mucus next to the epithelial surface. (Methods)
 
> [5.2] Bacteria were evaluated using fluorescence in situ hybridization, with bacterial 16S RNA probes and related to the periodic acid Schiff alcian blue stain. (Abstract)
 
`6.` Ringel Y, Maharshak N, Ringel-Kulka T, Sartor RB, Carroll I. High throughput sequencing reveals distinct microbial populations within the mucosal and luminal niches in healthy individuals. *Gut Microbes.* 2015;6:173 - 181. doi:10.1080/19490976.2015.1044711
 
> [6.1] Compared to mucosal samples, fecal samples had a signi fi cantly higher abundance of Firmicutes (41.4% vs. 29.1%, FDR < 0.0001, respectively) and of Actinobacteria (22% vs. 12.6%, FDR < 0.0001, respectively) compared to mucosal samples, and a signi fi cantly lower abundance of Bacteroidetes (20.2% vs. 26.3%, FDR < 0.05, respectively), and Proteobacteria (9.3% vs. 19.3%, FDR < 0.0001, respectively). (Discussion)
 
> [6.2] The dominant phyla in the LM and MAM were significantly different: Firmicutes (41.4% vs. 29.1%, FDR < 0.0001, respectively), Bacteroidetes (20.2% vs. 26.3%, FDR < 0.05, respectively), Actinobacteria (22% vs. 12.6%, FDR < 0.0001, respectively) and Proteobacteria (9.3% vs. 19.3%, FDR < 0.0001, respectively). (Abstract)
 
> [6.3] Methods: Fresh fecal samples and distal colonic mucosal biopsies collected from 24 healthy subjects before (fecal) and during (mucosa) a flexible sigmoidoscopy of an un-prepared bowel. (Abstract)
 
> [6.4] Results: LM and MAM populations were significantly different (ANOSIM: R = 0.49, P = 0.001). (Abstract)
 
`7.` Borgo F, Garbossa S, Riva A, et al. Body Mass Index and Sex Affect Diverse Microbial Niches within the Gut. *Frontiers in Microbiology.* 2018;9. doi:10.3389/fmicb.2018.00213
 
> [7.1] MAM sample analysis revealed enrichment in aerotolerant Proteobacteria, probably selected by a gradient of oxygen that decreases from tissue to lumen, and in Streptococcus and Clostridium spp. (Abstract)
 
> [7.2] In the present study, we compared LAM and MAM, by collecting stools and sigmoid brush samples of forty adults without large-bowel symptoms, and through a 16S rRNA gene next-generation sequencing (NGS) approach. (Abstract)
 
`8.` Hong P, Croix JA, Greenberg E, Gaskins H, Mackie R. Pyrosequencing-Based Analysis of the Mucosal Microbiota in Healthy Individuals Reveals Ubiquitous Bacterial Groups and Micro-Heterogeneity. *PLoS ONE.* 2011;6. doi:10.1371/journal.pone.0025042
 
> [8.1] The data demonstrate that mucosal-associated microbiota is comprised of Firmicutes (50.9% 621.3%), Bacteroidetes (40.2% 623.8%) and Proteobacteria (8.6% 64.7%), and that interindividual differences were apparent. (Body)
 
> [8.2] Phylogenetic affiliations Average percent abundance with respect to total microbial community± standard deviation across the three sampled sites in the same individual ABCDE Firmicutes Blautia 0.1660.05 0.16 60.05 0.50 60.16 0.53 60.29 0.27 60.24 Faecalibacterium 0.9960.85 1.31 60.39 0.38 60.15 0.89 60.47 0.43 60.54 Lactococcus 5.5562.33 1.40 60.77 3.37 61.71 3.24 62.98 6.35 63.83 Leuconostoc 15.61610.11 4.56 62.23 10.33 65.03 9.98 68.62 18.79 610.24 Streptococcus 0.6960.57 0.23 60.16 0.84 60.35 0.50 60.46 1.43 60.42 Weissella 20.64613.48 6.35 63.35 13.67 66.50 13.32 612.11 24.12 613.16 Veillonella 0.4160.38 0.14 60.12 0.38 60.11 0.32 60.33 0.57 60.32 Unclassified Clostridiales 3.76 62.57 4.42 60.83 6.96 60.95 6.95 64.25 2.85 62.56 Bacteroidetes Bacteroides 28.07620.75 26.93 613.94 41.18 611.60 18.39 69.37 15.26 621.65 Proteobacteria Acidovorax 0.0760.07 0.02 60.02 0.09 60.10 0.03 60.04 0.10 60.10 Acinetobacter 1.8061.24 0.54 60.35 1.50 60.78 1.29 61.46 3.37 62.23 Unclassified Enterobacteriaceae 3.52 60.87 0.68 60.47 3.37 61.25 1.39 61.32 3.81 61.64 Unclassified Pasteurellaceae 0.31 60.33 0.61 60.23 2.19 60.94 0.07 60.05 1.26 61.61 Ubiquity was defined based on the presence of related OTUs in (Results)
 
> [8.3] Sample collection F i v eh e a l t h ys u b j e c t s( 2w o m e n ,b o t ho fa g e5 0 ,a n d3m e n ,o fa g e s 50, 73 and 74), undergoing routine screening colonoscopy at Carle Foundation Hospital (Urbana, IL USA), were recruited for this study. (Methods)
 
> [8.4] Besides Firmicutes and Bacteroidetes, the abundance of Proteobacteria (1.7% to 19.2% of total microbial community) also varied inter-individually, but was generally much lower in abundance compared to Firmicutes and Bacteroidetes (Figure 1). (Results)
 
> [8.5] Differences in mucosal-associated microbiota along the colon Distinct differences in the taxonomical profiles of mucosal- associated microbiota and pooled stool microbiota resulted in two groups that clustered apart on the multidimensional scaling plot (MDS) regardless of the host origins (Figure 2), suggesting that the mucosal-associated microbiota examined in this study were not contaminated with feces. (Results)
 
`9.` Vuik F, Dicksved J, Lam S, et al. Composition of the mucosa-associated microbiota along the entire gastrointestinal tract of human individuals. *United European Gastroenterology Journal.* 2019;7:897 - 907. doi:10.1177/2050640619852255
 
> [9.1] The bacterial composition was characterised using 16 S rRNA sequencing with Illumina Miseq. (Abstract)
 
`10.` Lavelle A, Lennon G, Docherty N, et al. Depth-Dependent Differences in Community Structure of the Human Colonic Microbiota in Health. *PLoS ONE.* 2013;8. doi:10.1371/journal.pone.0078835
 
> [10.1] Mucosal biopsies sample the mucus gel layer, epithelium and variable amounts of the submucosa. (Results)
 
`11.` Miyauchi E, Taida T, Kawasumi M, Ohkusa T, Sato N, Ohno H. Analysis of colonic mucosa-associated microbiota using endoscopically collected lavage. *Scientific Reports.* 2022;12. doi:10.1038/s41598-022-05936-y
 
> [11.1] The values were calculated per ml for lavage and per sample for brush and biopsy. (Results)
 
> [11.2] A weighted UniFrac distance from the biopsy samples was much higher for feces compared to lavage and brush (Fig.2a), indicating that the composition of the mucosa-associated microbiota is distinct from that of the luminal microbiota, as shown in previous studies12,13. (Results)
 
> [11.3] The lavage samples showed an intermediate pattern in these data (Fig.2a–c). (Results)
 
> [11.4] Among the different types of colonic mucosal samples collected from healthy volunteers, the lavage samples contained a higher amount of bacterial DNA and were less contaminated with host DNA compared to mucosal brushing (brush) and biopsy. (Abstract)
 
`12.` Juge N. Relationship between mucosa-associated gut microbiota and human diseases. *Biochemical Society Transactions.* 2022;50:1225 - 1236. doi:10.1042/bst20201201
 
> [12.1] Given the limited number and size of studies and the heterogeneity of analysed samples, it is not possible to define MAM composition in healthy individuals; what is clear is that the composition of the MAM differs from that of the LM and the FM. (Results)
 
> [12.2] These inter-study discrepancies could be due to differences in the site and method of sampling as a recent study comparing microbial communities from aspirates, biopsies or stools showed a gradual change in the Bacteroidetes and Firmicutes ratio from biopsies to aspirates to stools [56]. (Results)
 
`13.` Swidsinski A, Weber J, Loening-Baucke V, Hale L, Lochs H. Spatial Organization and Composition of the Mucosal Flora in Patients with Inflammatory Bowel Disease. *Journal of Clinical Microbiology.* 2005;43:3380 - 3389. doi:10.1128/jcm.43.7.3380-3389.2005
 
> [13.1] of patients with concn of adherent bacteria greater than 10 10 /ml | 18 | 14 | 6 | 2 | 0 | 0 | | Mean concn of mucosal bacteria (10 10 /ml) ± SEM b , c | 9.1 ± 20.2 | 0.26 ± 0.31 | 0.31 ± 0.16 | 0.09 ± 0.18 | 0.02 ± 0.06 | 0.02 | | Ileum | 12.4 | 0.28 | 0.25 | 0.15 | 0.004 | | | Ascending colon | 6.8 | 0.21 | 0.39 | 0.04 | 0.02 | | | Sigmoid colon | 8.2 | 0.3 | 0.29 | 0.07 | 0.04 | | | Mean concn of mucosa adherent bacteria (10 10 /ml) c | 11.2 | 0.31 | 0.35 | 0.07 | 0.02 | 0.02 | | Mean concn of mucus spread bacteria (10 10 /ml) c | 2.4 | 0.08 | 0.34 | 0.11 | 0.03 | 0 | | Mean concn of mucus ceiling bacteria (10 10 /ml) c | 13.9 | 0.41 | 0.23 | 0.10 | 0.02 | > 0 a | | Mean % of intact epithelial surface covered with mucosal bacteria c | 85 | 52 | 52 | 21 | 7 | Unclear | (Results)
 
> [13.2] Mucosal bacteria were found at concentrations greater than 10(9)/ml in 90 to 95% of IBD patients, 95% of patients with self-limiting colitis, 65% of IBS patients, and 35% of healthy controls. (Abstract)
 
`14.` Shen XJ, Rawls J, Randall T, et al. Molecular characterization of mucosal adherent bacteria and associations with colorectal adenomas. *Gut microbes.* 2010;1:138 - 147. doi:10.4161/gmic.1.3.12360
 
> [14.1] Overall, Firmicutes (62%), Bacteroidetes (26%) and Proteobacteria (11%) were the most dominant phyla. (Abstract)
 
`15.` Mukhopadhya I, Martin JC, Shaw S, Mckinley A, Gratz S, Scott K. Comparison of microbial signatures between paired faecal and rectal biopsy samples from healthy volunteers using next-generation sequencing and culturomics. *Microbiome.* 2022;10. doi:10.1186/s40168-022-01354-4
 
> [15.1] However, the rectal biopsy tissue (B) samples had a significantly altered bacterial signature with greater abundance of Proteobacteria and Acidobacteria compared to faecal (F) and faecal homogenised (FHg) samples. (Results)
 
`16.` Yilmaz B, Moulin S, Heimgartner B, et al. Spatial mapping of human colonic niches reveals rapid, mucus-specific microbiota disruption after bowel cleansing. *Gut Microbes.* 2026;18. doi:10.1080/19490976.2026.2635866
 
> [16.1] In contrast, the aspirated mucus samples captured a distinct community enriched for Enterobacteriaceae, suggesting that facultative anaerobes preferentially colonize the outer mucus layer rather than the deeper mucosal surface. (Body)
 
> [16.2] Critically, spatially resolved sampling revealed a key refinement: the Enterobacteriaceae expansion was confined almost exclusively to the superficial mucus layer, a glycan-rich, dynamically oxygenated compartment that is particularly susceptible to mechanical disturbance during lavage, whereas deeper mucus and mucosa-associated communities remained comparatively stable. (Abstract)
 
`17.` Flynn KJ, Ruffin M, Turgeon D, Schloss P. Spatial Variation of the Native Colon Microbiota in Healthy Adults. *Cancer Prevention Research.* 2017;11:393 - 402. doi:10.1158/1940-6207.capr-17-0370
 
> [17.1] The distal mucosa had higher relative abundance of Finegoldia, Murdochiella, Peptoniphilus, Porphyromonas, and Anaerococcus The proximal mucosa had more of the genera Enterobacteriaceae, Bacteroides, and Pseudomonas The classification model performed well when classifying mucosal samples into proximal or distal sides (AUC = 0.808). (Abstract)
 
`18.` De Cárcer DA, Cuív PÓ, Wang T, et al. Numerical ecology validates a biogeographical distribution and gender-based effect on mucosa-associated bacteria along the human colon. *The ISME Journal.* 2011;5:801-809. doi:10.1038/ismej.2010.177
 
> [18.1] Conversely, the analyses suggest that several members of the Enterobacteriaceae increase in relative abundance towards the rectum. (Abstract)
 
`19.` Mylonaki M, Rayment N, Rampton D, Hudspith B, Brostoff J. Molecular Characterization of Rectal Mucosa‐Associated Bacterial Flora in Inflammatory Bowel Disease. *Inflammatory Bowel Diseases.* 2005;11:481–487. doi:10.1097/01.mib.0000159663.62651.4f
 
> [19.1] Epithelium-associated counts of bifidobacteria in active [median 15/mm of epithelial surface (range, 4-56), n = 14] and quiescent ulcerative colitis [26/mm (range, 11-140), n = 19] were lower than in controls [56/mm (range, 0-144), n = 14; P = 0.006 and P = 0.03, respectively]. (Abstract)
 
`20.` Matsumoto H, Kuroki Y, Higashi S, et al. Analysis of the colonic mucosa associated microbiota (MAM) using brushing samples during colonic endoscopic procedures.. *Journal of clinical biochemistry and nutrition.* 2019;65 2:132-137. doi:10.3164/jcbn.19-3
 
> [20.1] The aim was to investigate mucosa/hyphenminusassociated microbiota using brush samples during endoscopic procedures and compare with fecal microbiota. (Results)
 
`21.` Lavelle A, Lennon G, O'Sullivan Ó, et al. Spatial variation of the colonic microbiota in patients with ulcerative colitis and control volunteers. *Gut.* 2015;64:1553 - 1561. doi:10.1136/gutjnl-2014-307873
 
`22.` Lennon G, Bambury N, Lavelle A, et al. Correlations between colonic crypt mucin chemotype, inflammatory grade and Desulfovibrio species in ulcerative colitis. *Colorectal Disease.* 2014;16. doi:10.1111/codi.12503
 
`23.` Swidsinski A, Ladhoff A, Pernthaler A, et al. Mucosal flora in inflammatory bowel disease.. *Gastroenterology.* 2002;122 1:44-54. doi:10.1053/gast.2002.30294
 
`24.` Ahmed S, Macfarlane G, Fite A, McBain A, Gilbert P, Macfarlane S. Mucosa-Associated Bacterial Diversity in Relation to Human Terminal Ileum and Colonic Biopsy Samples. *Applied and Environmental Microbiology.* 2007;73:7435 - 7442. doi:10.1128/aem.01143-07
 
`25.` Macfarlane S, Woodmansey EJ, Macfarlane G. Colonization of Mucin by Human Intestinal Bacteria and Establishment of Biofilm Communities in a Two-Stage Continuous Culture System. *Applied and Environmental Microbiology.* 2005;71:7483 - 7492. doi:10.1128/aem.71.11.7483-7492.2005
 
`26.` Pullan R, Thomas G, Rhodes M, et al. Thickness of adherent mucus gel on colonic mucosa in humans and its relevance to colitis.. *Gut.* 1994;35:353 - 359. doi:10.1136/gut.35.3.353
 
`27.` Martinez-Medina JN, Ghazisaeedi F, Kramer C, et al. Mucosal washes are useful for sampling intestinal mucus-associated microbiota despite low biomass. *Gut Microbes.* 2023;17. doi:10.1080/19490976.2025.2464296
