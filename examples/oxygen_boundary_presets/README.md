# Epithelial Oxygen Boundary & Anaerobic Maintenance Presets (`GAP-08` / Issue `#470`)

Shipped `GutIBM` defaults preserve `oxygen.enabled = false` and `oxygen.epithelial_conc = 55.0e-6 mol/m^3` (`~0.042 mmHg`) for backward-compatible baseline trajectories. As established in `docs/CYCLE2_EMPIRICAL_CALIBRATION_REVIEW.md` (§3, `GAP-08`), choosing an epithelial oxygen value requires specifying **which physical plane** `z = 0` represents rather than applying a blind `1,000x` multiplier:

1. **Apical mucus-base Dirichlet concentration (`oxygen.epithelial_boundary = "dirichlet"`):** `oxygen.epithelial_conc` prescribes dissolved `O2` *after* colonocyte mitochondrial consumption right at the apical epithelial surface (`< 1 mmHg` in homeostatic human colon; `~1–5 mmHg` in murine colonic surface/crypt micro-niches).
2. **Submucosal tissue-reservoir Robin delivery (`oxygen.epithelial_boundary = "robin"`):** `oxygen.epithelial_conc` prescribes the vascularized lamina propria reservoir (`~38–42 mmHg`, `5.0e-2–5.5e-2 mol/m^3`) coupled across finite epithelial permeability `oxygen.epithelial_transfer_coeff` (`k = 1.2e-6 m/s` homeostatic; elevated under inflammation), with `oxygen.z_gradient = false` required by `InputParser`.

## Preset Matrix

| File | Boundary Mode (`z = 0`) | `oxygen.epithelial_conc` (`mol/m^3`) | `transfer_coeff` (`m/s`) | `vbf_sink` (`1/s`) | Host / Compartment Prior |
| :--- | :--- | ---: | ---: | ---: | :--- |
| `apical_dirichlet_human_anoxic.json` | `dirichlet` | `5.5e-5` (`0.042 mmHg`) | `0.0` | `1.0e-3` | Homeostatic human near-anoxic apical mucus base |
| `apical_dirichlet_murine_surface.json` | `dirichlet` | `3.93e-3` (`3.0 mmHg`) | `0.0` | `1.0e-3` | Murine distal colonic surface / crypt prior (Konjar 2021; Albenberg 2014) |
| `tissue_robin_healthy.json` | `robin` | `5.5e-2` (`42 mmHg`) | `1.2e-6` | `1.0e-3` | Submucosal tissue reservoir + intact colonocyte barrier (`funded` delivery) |
| `tissue_robin_inflamed.json` | `robin` | `5.5e-2` (`42 mmHg`) | `3.6e-6` | `2.0e-4` | Hyperemic / leaky inflamed mucosa + depleted anaerobe scavenging |

## Anaerobic Carbon Maintenance Priors (`oxygen.anaerobic_maintenance_factor`)

- **Shipped ceiling (`15.0`):** Derived from idealized aerobic-to-fermentative ATP yield stoichiometry (`~30 / 2` mol ATP per mol glucose; inert when `carbon_maintenance_rate = 0.0`).
- **Empirical / flux-analysis sensitivity range (`3.4x–8.0x`, midpoint `5.5x` used here):** Chemostat and metabolic flux comparisons (`Boecker et al. 2021`, `Carlson et al. 2026` conditional scenarios) place realized anaerobic-to-aerobic specific glucose maintenance ratios between `3.4x` and `8.0x` (`carbon_maintenance_rate = 2.1e-5 mol/(s*kg wet biomass)` aerobic anchor).
