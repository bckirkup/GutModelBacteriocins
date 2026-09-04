---
name: searching-literature-evidence
description: Search the peer-reviewed literature with the Consensus MCP server to source a GutIBM constant — mucus rheology, diffusion and binding coefficients, mucosal densities, degradation rates — including how a hit becomes a provenance record in types.h and docs/PARAMETERS.md. Use whenever a physical or biological constant needs a citation, or when asked what the literature says about a mechanism. Pairs with the org-level consensus-literature-retrieval skill, which owns retrieval mechanics.
---

# Searching the Literature (Consensus MCP)

## Retrieval mechanics are in the org-level skill

Load `consensus-literature-retrieval` (`~/.agents/skills/`) before searching. It
owns the tool surface, `include_full_text_chunks: true` — which is mandatory and
returns Results, Methods and tables, including for paywalled articles — query
construction, filter behaviour, result handling, and recording which section of
the paper a number was read from.

This skill is the other half: what needs sourcing in GutIBM, and what a hit is
allowed to become here.

## Query construction

A useful query names the measured quantity plus its preparation:

- Good: `mucin gel diffusion coefficient FRAP protein retardation`
- Weak: `how fast do colicins move through mucus`

Quantities this repo actually needs sourced, and the words that find them:

- Mucus rheology — `yield stress`, `viscoelastic`, `oscillatory rheometry`,
  `cecal`, `colonic mucus`.
- Mucosal *Enterobacteriaceae* abundance — `CFU per gram`, `mucosa-associated`,
  `qPCR 16S copies`, `biopsy`.
- Bacteriocin transport — `diffusion coefficient`, `mucin binding`,
  `hindered diffusion`, `colicin`.
- Receptor kinetics — `TonB-dependent transporter`, `Kd`, `binding affinity`,
  `expression level`.
- Degradation — `protease`, `half-life`, `stability`, `inactivation rate`.

The paper establishing that mucin retards a solute is rarely the one that
measured the coefficient.

## Filter discipline

Two filters are actively harmful here:

- `medical_mode=true` restricts to ~8M top medical documents and drops
  *Applied and Environmental Microbiology*, *Journal of Bacteriology*,
  *Biophysical Journal*, and the rheology literature — that is, nearly every
  journal a GutIBM constant comes from.
- `human=true` drops in-vitro, murine and gnotobiotic work, which is where
  essentially all mucus-layer and colicin measurements were made.

Reasonable filters: `domain="bio,med,chem"` when a query drags in unrelated
fields. Do not set `year_min` — a 1990s diffusion measurement is not stale.

## Capture the unit and the assay, not just the number

This repository has a units audit (`docs/UNITS_AUDIT.md`) because unit and
assay confusion is the recurring defect, not arithmetic. A density reported as
CFU/g wet weight, copies/g dry weight, or cells/mL of luminal contents are three
different numbers under one heading. Record which one the paper measured, in the
paper's own units, and convert visibly in the code.

Watch the definition, not just the magnitude: total versus mucosa-associated
abundance, or bulk versus near-wall viscosity, differ by orders of magnitude
under similar titles. Taking a headline figure without checking what was
measured is the same error as having no source.

## From a hit to a sourced constant

Provenance in this repo currently lives mostly in `docs/PARAMETERS.md` — for
example *"Yield stress of mucus | 46 ± 9 Pa | Greter et al. 2026, cecal
rheology"* — while the definitions in `src/core/types.h` carry only a unit:

```cpp
static constexpr Real GUT_PH = 6.8;   // colonic pH
```

When you source a constant, put the citation **at the definition** as well as in
the table, because the table drifts and the header is what the next reader sees:

```cpp
// Proximal colon luminal pH, in vivo radiotelemetry, mean 6.8 (range 5.7-7.4).
// <Author> et al. <year>, <journal> (DOI: <doi>). Grade A.
static constexpr Real GUT_PH = 6.8;
```

State **what was measured**, **in what preparation**, the value with its range,
and author + year + journal + DOI. Then grade it:

- **A** — direct measurement of this quantity in this setting (human colon,
  mucus layer, the organism modelled).
- **B** — direct measurement in an analogous setting (murine cecal mucus for
  human colonic mucus; a related bacteriocin for colicin).
- **C** — inferred, estimated, or a declared assumption. Every Grade C constant
  is a standing liability; list it in `docs/PARAMETERS.md` as assumed.

If no source exists, say so explicitly rather than inventing a plausible number.
`k_ROS` and `b12_initial_conc` are known-uncited values: a declared Grade C with
a sensitivity sweep is honest, a fabricated citation is not.

## What this search must never be used for

Do not search for a value that makes a calibration or validation target come out
right — the density reconciliation, the Section 8 validation, or an AWS
calibration run. Sourcing a constant independently is what makes the later
comparison a real test; screening candidate papers by which value helps destroys
that test just as surely as fitting the constant by hand.

If a sourced constant makes a target worse, that is a result: report it.
