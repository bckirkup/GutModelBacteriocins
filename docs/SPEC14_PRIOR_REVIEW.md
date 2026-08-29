# Spec 14 priors: lysogen prevalence, and the induction rate that bounds it

Review of the Spec 14 rev3 prior revision (lysogen prevalence re-sourced from
Kim & Bae 2018 to Dikareva et al. 2023, default 0.8, sweep {0.3, 0.8, 1.0}, and
the restored genomic-vs-inducible caveat), plus one independent source that the
revision does not yet account for.

Nothing here unblocks implementation. `docs/SPEC14_SECTION8_VALIDATION.md`
found that GutIBM does not reproduce the in-vitro exclusion Spec 14 exists to
explain, and that failure is upstream of every parameter discussed below. This
document exists so that the prior work is not lost, and because one of the
constraints found here would change the campaign design even if §8 were passing.

## 1. The re-sourcing holds, with a cohort caveat

Checked against the primary sources rather than the review's paraphrase.

Dikareva et al. 2023 (*Front Microbiol* 14:1254535) is the right source: human
fecal, 6,186 MAGs, 7,165 prophage sequences >10 kb, **>70% of near-complete MAGs
are lysogens**, and Enterobacteriaceae is named explicitly among the
highest-prevalence families alongside Oscillospiraceae and Enterococcaceae.
Right host, right clade, and three orders of magnitude more genomes than the 4
bins behind the Kim & Bae figure. A default of 0.8 for our clade is better
supported after the swap than before it.

One caveat the revision does not carry: the cohort is **infant-weighted** — 88
infants at 3 weeks / 3 / 6 / 12 months (323 samples) against their parents (138
samples) — and the paper itself reports that for several families prophage
prevalence is *higher* at the early infant timepoints than in adults. The
family-level Enterobacteriaceae statement is not broken out by age in the
abstract. This is a second-order concern next to §2 below, but it means 0.8 is
an upper-ish estimate for the adult gut, not a central one.

## 2. Prevalence is not independently free: the observable is the product

This is the finding that would change the campaign.

The Current Biology 2025 paper flagged for later reading —
*Abundance measurements reveal the balance between lysis and lysogeny in the
human gut microbiome*, doi:10.1016/j.cub.2025.03.073 — estimates that **phage
induction and lysis occur at ~0.001–0.01 events per bacterium per day** in the
human gut, from phage-particle-to-cell (~1:100) and phage-genome-to-bacterial-
genome (~4:1) ratios combined with imaging. That is a *community-average
realized* rate, so it measures

```text
observed rate  =  P(inducible lysogen)  x  induction rate per lysogen
```

— the **product** of the two quantities Spec 14 treats as separate parameters.

Spec 14's `spontaneous_induction_rate = 1.75e-7 /s` is 1.5% per day per lysogen
(that is where the number came from: `-ln(1-0.015)/86400`). Multiplying through
the prevalence sweep:

| `lysogen_prevalence_init` | implied community induction, per bacterium per day | vs. observed 0.001–0.01 |
|---|---|---|
| 0.3 | 0.0045 | in band, upper half |
| **0.8 (default)** | 0.0121 | **above the ceiling** |
| 1.0 | 0.0151 | **1.5× the ceiling** |

So the sweep as specified does not hold the observable fixed while varying a
free parameter — it **walks the model's induction rate from the top of the
observed band to 1.5× above it**, and it does so before Change 1's damage-driven
Hill term adds anything. The observed band covers *total* induction and lysis,
not the spontaneous component alone, so the damage term has to fit inside the
same ceiling that the spontaneous baseline is already at.

Two consequences for the design:

1. The prevalence sweep should be run with `spontaneous_induction_rate`
   compensating, so that the product stays inside 0.001–0.01/bacterium/day.
   Otherwise "prevalence" and "induction rate" are one axis sampled twice, and
   any prevalence effect the campaign reports is confounded with a total-lysis
   effect that is independently constrained.
2. Realized total induction per bacterium per day becomes a **validation
   output**, reported alongside the campaign result and checked against the
   band — the same discipline `docs/SPEC13_LYSIS_SELECTION.md` already applies
   to realized `mortality_lysis / divisions`.

Strength of the constraint: this is a soft prior, not a hard gate. The estimate
is a community average over a Bacteroidales-dominated microbiome, and
Enterobacteriaceae λ lysogens in vivo may legitimately sit above it — De Paepe's
central finding is precisely that reactivation is frequent in the mouse gut. But
it is the **independent replication our parameterization otherwise lacks**: it
does not come from Petit's group, and it is the first constraint on our phage
numbers that is not correlated with the two sources supplying burst size,
lysogenization fraction, and spontaneous rate.

The same paper also reports that most gut phage are effectively temperate and
that free particles run ~1:100 against cells, which is independent support for
keeping Change 3 (the free-phage field) deferred: at ~1 particle per 100 cells,
a continuum field is representing a quantity that is nearly always zero in any
voxel we resolve.

## 3. `is_lysogen` means inducible, and the asymmetry cuts against the mechanism

The restored caveat is correct and load-bearing: every prevalence number above
comes from **detecting integrated prophage sequences in genomes**, and many
detected prophages are cryptic or defective remnants that cannot excise.
Genomic carriage is an upper bound on inducible carriage. `is_lysogen` in the
model must mean *inducible* lysogen, and 0.3 is therefore a plausible
"most detected prophages are cryptic" scenario rather than a low-end strawman.

What the revision does not draw out is that the two arms of the rev2 argument
do not degrade together. Rev2 argued that high prevalence makes most competitors
inducible **and** most targets superinfection-immune. Under the cryptic reading
only the first of those weakens: **a cryptic prophage still confers repressor-
mediated immunity to its own phage even though it cannot excise.** So as the
inducible fraction falls from 0.8 to 0.3, the kin-lysogenization channel and the
second killing channel weaken while the immunity blockade stays roughly where
genomic carriage puts it. The mechanism is therefore *less* consequential at the
low end of the sweep than a symmetric reading of the prevalence number suggests,
and the sweep's endpoints are not two settings of one dial — they are two
different structural claims about the community.

If Change 2 is ever implemented, this argues for representing immunity and
inducibility as **separate flags** rather than deriving both from `is_lysogen`.
That costs one bit per agent and nothing else, and it is much cheaper to do at
implementation time than to retrofit.

## 4. Sequence-detectable is not functionally active

This is now the third instance of one failure mode, and it should be tracked as
a class rather than rediscovered:

| Instance | Molecular measure | Functional quantity | Ratio |
|---|---|---|---|
| Mucosal density (`docs/gutibm_density_reconciliation.md`) | 16S qPCR copies/mg | culturable CFU | 100–1000× |
| Lysogen prevalence (this document) | integrated prophage in MAGs | inducible prophage | unquantified; strictly ≤ 1 |
| Free phage vs. induction (§2) | phage genomes per bacterial genome, ~4:1 | induction events/bacterium/day, ~0.001–0.01 | the genome ratio does not predict the rate |

In each case a sequencing-derived number was available and precise, the
functional number was not, and the sequencing number was substituted. The
generalizable rule: **when a parameter is defined functionally (inducible,
viable, active) but sourced from a sequence census, record it as an upper bound
and put the functional fraction in the sweep** — do not adopt the census value
as the central estimate. Added to `AGENTS.md`.

## 5. Provenance status

The correlated-provenance caveat improves but does not clear. De Paepe 2016
(burst size 12.1, lysogenization fraction 0.19, spontaneous rate) and Henrot
2026 (bacteriocin specificity, the in-vitro/in-vivo reversal) still share a
senior author, and they remain the sources for every *mechanism* parameter.
Dikareva 2023 and Current Biology 2025 are independent, and they constrain
prevalence and the induction rate respectively — i.e. the two quantities that
scale the mechanism, not the mechanism itself. That is a real improvement in the
part of the parameterization most likely to be load-bearing, and it should be
stated that way rather than as "the provenance concern is resolved".
