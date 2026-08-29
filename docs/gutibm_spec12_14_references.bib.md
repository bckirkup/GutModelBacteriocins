# References — Spec 12 Addendum A and Spec 14

BibTeX from the two literature reviews backing these specs.
Sources: task 1114ae84-1eb2-4f68-88f3-378a90ee6d96 (Damkohler anchors),
task 44e9e0b4-233f-4123-8301-93399cb96056 (prophage induction).

## Primary motivating preprints

```bibtex
@article{henrot2026dnadamaging,
  author={Henrot, Caroline and Debarbieux, Laurent and Petit, Marie-Agn\`es},
  title={DNA-damaging bacteriocins from human Escherichia coli intestinal isolates
         trigger prophage induction and promote lysogeny},
  journal={bioRxiv}, year={2026}, month={May},
  doi={10.64898/2026.05.26.727859},
  note={Preprint. Motivates Spec 14.}
}
@article{scheidweiler2026oxygen,
  author={Scheidweiler, David and Bornet, Elise and Ochner, Hannah and Blasche, Sonja
          and Thiriet-Rupert, Stanislas and Dorison, Louis and Beloin, Christophe
          and Davit, Yohan and Gobaa, Samy and Bharat, Tanmay A. M.
          and Westermann, Alexander J. and Patil, Kiran R. and Ghigo, Jean-Marc},
  title={Oxygen gradients reshape cross-feeding through emergent spatial organization
         of gut commensal bacteria},
  journal={bioRxiv}, year={2026},
  doi={10.64898/2026.05.02.721930},
  note={Preprint. Supports the Damkohler mechanism QUALITATIVELY only -- reports no
        local packing density, calibrated per-cell uptake, or clonal aggregate radius.
        Do not cite as quantitative confirmation of our crossover.}
}
```

## Phage and prophage biology (Spec 14)

```bibtex
@article{paepe2016carriage,
  author={De Paepe, Marianne and Tournier, Laurent and Moncaut, Elisabeth and Son, Olivier
          and Langella, Philippe and Petit, Marie-Agn\`es},
  title={Carriage of $\lambda$ Latent Virus Is Costly for Its Bacterial Host due to
         Frequent Reactivation in Monoxenic Mouse Intestine},
  journal={PLOS Genetics}, year={2016}, volume={12}, pages={e1005861},
  doi={10.1371/journal.pgen.1005861},
  note={Source of the IN VIVO burst size (12.1 PFU/cell), lysogenization fraction
        (0.19), and spontaneous induction (1-2 percent of population). Use these,
        not in vitro values.}
}
@article{kourilsky1973lysogenization,
  author={Kourilsky, Philippe},
  title={Lysogenization by bacteriophage lambda},
  journal={Molecular and General Genetics}, year={1973}, volume={122}, pages={183--195},
  doi={10.1007/bf00435190}
}
@article{shao2008adsorption,
  author={Shao, Yongping and Wang, Ing-Nang},
  title={Bacteriophage Adsorption Rate and Optimal Lysis Time},
  journal={Genetics}, year={2008}, volume={180}, pages={471--482},
  doi={10.1534/genetics.108.090100}
}
@article{ellis1939growth,
  author={Ellis, Emory L. and Delbr\"uck, Max},
  title={The Growth of Bacteriophage},
  journal={Journal of General Physiology}, year={1939}, volume={22}, pages={365--384},
  doi={10.1085/jgp.22.3.365}
}
@article{nabergoj2018growthrate,
  author={Nabergoj, Dominik and Modic, Petra and Podgornik, Ale\v{s}},
  title={Effect of bacterial growth rate on bacteriophage population growth rate},
  journal={MicrobiologyOpen}, year={2018}, volume={7},
  doi={10.1002/mbo3.558}
}
@article{dikareva2023prophages,
  author={Dikareva, Evgenia and others},
  title={An extended catalog of integrated prophages in the infant and adult fecal
         microbiome shows high prevalence of lysogeny},
  journal={Frontiers in Microbiology}, year={2023}, volume={14}, pages={1254535},
  doi={10.3389/fmicb.2023.1254535},
  note={Human fecal, 6186 MAGs, 7165 prophage sequences; >70\% of near-complete MAGs
        are lysogens, Enterobacteriaceae among the highest-prevalence families.
        Supersedes the murine 4-bin Kim \& Bae 2018 figure as the prevalence anchor.
        Reports GENOMIC prophage carriage: an upper bound on inducible carriage, and
        the cohort is infant-weighted. See docs/SPEC14_PRIOR_REVIEW.md.}
}
@article{currentbiology2025lysislysogeny,
  title={Abundance measurements reveal the balance between lysis and lysogeny in the
         human gut microbiome},
  journal={Current Biology}, year={2025},
  doi={10.1016/j.cub.2025.03.073},
  note={Preprint doi:10.1101/2024.09.27.614587. Phage particles ~1:100 to cells,
        phage genomes ~4:1 to bacterial genomes, induction and lysis ~0.001--0.01 per
        bacterium per day. Independent of the Petit-group sources: constrains the
        PRODUCT of lysogen prevalence and per-lysogen induction rate, which Spec 14
        treats as two free parameters. See docs/SPEC14_PRIOR_REVIEW.md.}
}
```

## Oxygen physiology and transport (Addendum A)

```bibtex
@article{alexeeva2002oxygen,
  author={Alexeeva, Svetlana and Hellingwerf, Klaas J. and Teixeira de Mattos, M. Joost},
  title={Quantitative Assessment of Oxygen Availability: Perceived Aerobiosis and Its
         Effect on Flux Distribution in the Respiratory Chain of Escherichia coli},
  journal={Journal of Bacteriology}, year={2002}, volume={184}, pages={1402--1406},
  doi={10.1128/JB.184.5.1402-1406.2002}
}
@article{rice1978oxygenlimited,
  author={Rice, C. W. and Hempfling, W. P.},
  title={Oxygen-limited continuous culture and respiratory energy conservation in
         Escherichia coli},
  journal={Journal of Bacteriology}, year={1978}, volume={134}, pages={115--124},
  doi={10.1128/JB.134.1.115-124.1978}
}
@article{friedman2018microbes,
  author={Friedman, Elliot S. and Bittinger, Kyle and Esipova, Tatiana V. and others},
  title={Microbes vs. chemistry in the origin of the anaerobic gut lumen},
  journal={PNAS}, year={2018}, volume={115}, pages={4170--4175},
  doi={10.1073/pnas.1718635115}
}
@article{zheng2015hypoxia,
  author={Zheng, Leon and Kelly, Caleb J. and Colgan, Sean P.},
  title={Physiologic hypoxia and oxygen homeostasis in the healthy intestine},
  journal={American Journal of Physiology-Cell Physiology},
  year={2015}, volume={309}, pages={C350--C360},
  doi={10.1152/ajpcell.00191.2015}
}
@article{vasilakou2020ecolimetabolism,
  author={Vasilakou, Eleni and van Loosdrecht, Mark C. M. and Wahl, S. Aljoscha},
  title={Escherichia coli metabolism under short-term repetitive substrate dynamics:
         adaptation and trade-offs},
  journal={Microbial Cell Factories}, year={2020}, volume={19},
  doi={10.1186/s12934-020-01379-0}
}
```

## Spatial organization and packing density (Addendum A)

```bibtex
@article{swidsinski2005spatial,
  author={Swidsinski, Alexander and Weber, Jutta and Loening-Baucke, Vera
          and Hale, Laura P. and Lochs, Herbert},
  title={Spatial Organization and Composition of the Mucosal Flora in Patients with
         Inflammatory Bowel Disease},
  journal={Journal of Clinical Microbiology}, year={2005}, volume={43}, pages={3380--3389},
  doi={10.1128/JCM.43.7.3380-3389.2005}
}
@article{earle2015quantitative,
  author={Earle, Kristen A. and Billings, Gabriel and Sigal, Michael and others},
  title={Quantitative Imaging of Gut Microbiota Spatial Organization},
  journal={Cell Host and Microbe}, year={2015}, volume={18}, pages={478--488},
  doi={10.1016/j.chom.2015.09.002}
}
@article{welch2017spatial,
  author={Mark Welch, Jessica L. and Hasegawa, Yuko and McNulty, Nathan P.
          and Gordon, Jeffrey I. and Borisy, Gary G.},
  title={Spatial organization of a model 15-member human gut microbiota established
         in gnotobiotic mice},
  journal={PNAS}, year={2017}, volume={114}, pages={E9105--E9114},
  doi={10.1073/pnas.1711596114}
}
@article{mondragon2022threedimensional,
  author={Mondrag\'on-Palomino, Octavio and Poceviciute, Roberta and Lignell, Antti and others},
  title={Three-dimensional imaging for the quantification of spatial patterns in
         microbiota of the intestinal mucosa},
  journal={PNAS}, year={2022}, volume={119},
  doi={10.1073/pnas.2118483119}
}
```

## Reaction-diffusion and biofilm oxygen limitation (Addendum A)

```bibtex
@article{wessel2014oxygen,
  author={Wessel, Aimee K. and Arshad, Talha A. and Fitzpatrick, Mignon and
          Connell, Jodi L. and Bonnecaze, Roger T. and Shear, Jason B. and Whiteley, Marvin},
  title={Oxygen Limitation within a Bacterial Aggregate},
  journal={mBio}, year={2014}, volume={5},
  doi={10.1128/mBio.00992-14},
  note={Closest direct comparison: dense P. aeruginosa aggregate at ~1e12 cells/mL,
        predicted depletion radius ~35 um, O2 penetration ~60 um.}
}
@article{stewart2016reactiondiffusion,
  author={Stewart, Philip S. and Zhang, Tianyu and Xu, Ruifang and others},
  title={Reaction-diffusion theory explains hypoxia and heterogeneous growth within
         microbial biofilms associated with chronic infections},
  journal={npj Biofilms and Microbiomes}, year={2016}, volume={2},
  doi={10.1038/npjbiofilms.2016.12}
}
@article{wu2018hypoxia,
  author={Wu, Yilin and Klapper, Isaac and Stewart, Philip S.},
  title={Hypoxia arising from concerted oxygen consumption by neutrophils and
         microorganisms in biofilms},
  journal={Pathogens and Disease}, year={2018}, volume={76}, pages={fty043},
  doi={10.1093/femspd/fty043}
}
@article{vandenberg2021diffusion,
  author={van den Berg, Lenno and van Loosdrecht, Mark C. M. and de Kreuk, Merle K.},
  title={How to measure diffusion coefficients in biofilms: A critical analysis},
  journal={Biotechnology and Bioengineering}, year={2021}, volume={118}, pages={1273--1285},
  doi={10.1002/bit.27650}
}
@article{lafitte2007diffusion,
  author={Lafitte, G\'eraldine and Thuresson, Krister and S\"oderman, Olle},
  title={Diffusion of nutrient molecules and model drug carriers through mucin layer
         investigated by magnetic resonance imaging with chemical shift resolution},
  journal={Journal of Pharmaceutical Sciences}, year={2007}, volume={96}, pages={258--263},
  doi={10.1002/jps.20749},
  note={Basis for treating glucose as only mildly retarded in mucus.}
}
```
