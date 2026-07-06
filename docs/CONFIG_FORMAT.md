# GutIBM Input Config Format

GutIBM simulation configs are **strict JSON objects**. They validate with standard
JSON tools (`python -m json.tool`, `jq`, IDE parsers) and are loaded by
`InputParser::parse()` via the built-in JSON document parser.

For **batch parameter scans** (multiple runs from one manifest), see
[BATCH_RUNNER.md](BATCH_RUNNER.md) — a separate JSON schema from simulation input.

## Comment convention

JSON does not support `//` or `/* */` comments. Use a metadata field instead:

```json
{
  "_comment": [
    "Human-readable notes about this scenario.",
    "Second line is fine — string or array of strings."
  ],
  "total_time": 86400
}
```

Keys starting with `_` are ignored by the parser (except as documentation).

## Top-level structure

| Key | Type | Description |
|-----|------|-------------|
| `_comment` | string or string[] | Documentation only (ignored) |
| `total_time`, `bio_dt`, … | number | Scalar parameters — see [PARAMETERS.md](PARAMETERS.md) |
| `initial_strains` | array | Strain objects with `type`, `count`, `plasmids`, … |
| `fixes` | string[] | Fix plugin names in execution order (empty = all registered) |
| `oxygen.enabled`, `acetate.enabled`, … | bool | Nested toggles via dot keys |
| `checkpoint_file` | string | HDF5 checkpoint for restart |
| `checkpoint_step` | string | Optional step group name to restore |

## Fix selection

Registered Fix plugins (default execution order when `fixes` is omitted):

```
metabolism → bacteriocin → receptor → motility → conjugation → cdi → mutation → mechanics
```

Override with a `fixes` array:

```json
{
  "fixes": ["metabolism", "bacteriocin", "receptor", "conjugation", "mutation"]
}
```

## Feature toggles (Spec 1 & 3)

Many subsystems use **dot-key toggles** in flat JSON:

| Key | Default | Subsystem |
|-----|---------|-----------|
| `oxygen.enabled` | false | Epithelial O₂ gradient, aerobic boost, ROS→SOS |
| `acetate.enabled` | false | Dynamic acetate production/scavenging |
| `mucin.enabled` | false | Mucin glycoprotein field + liberation |
| `protease.enabled` | true | Toxin protease decay in QSSA |
| `fur.enabled` | false | Fur-regulated iron receptor expression |
| `cdi.enabled` | false | Contact-dependent inhibition |
| `motility.enabled` | false | Run-and-reverse swimming |
| `crypts_enabled` | false | Crypt refugia (zero-flow zones) |
| `adaptive_dt_enabled` | false | CFL-like adaptive biological timestep |
| `use_fmm` | false | Barnes–Hut FMM far-field acceleration |
| `gpu_enabled` | false | CUDA GPU path (CUDA build required) |

Full parameter lists: [PARAMETERS.md](PARAMETERS.md).

## Per-strain fields

`initial_strains[]` objects support simulation-specific overrides:

```json
{
  "type": 1,
  "count": 120,
  "mu_max": 5.5e-4,
  "plasmids": ["ColE1", "ColB"],
  "conjugative": true,
  "cdi_type": 1,
  "cdi_immunity": 1
}
```

## Example (minimal)

```json
{
  "_comment": "Resident colony with colicin production",

  "total_time": 86400,
  "bio_dt": 60,
  "seed": 42,

  "domain_x": 0.001,
  "domain_y": 0.001,
  "domain_z": 0.0001,
  "grid_dx": 2e-6,

  "initial_strains": [
    {
      "type": 1,
      "count": 500,
      "mu_max": 5.5e-4,
      "plasmids": ["ColE1", "ColB"],
      "conjugative": true
    }
  ],

  "hdf5_file": "output.h5",
  "hdf5": {
    "schedule": {
      "summary": 60,
      "agents": 60,
      "grid": 0,
      "lineage": 100,
      "genome": 100
    }
  }
}
```

## Checkpoint restart

```json
{
  "checkpoint_file": "checkpoint.h5",
  "checkpoint_step": "step_000100",
  "total_time": 172800
}
```

## Validation

From the repo root:

```bash
bash scripts/validate_config_json.sh
```

CI runs this on every push to ensure examples and parser fixtures remain valid JSON.

## Legacy line-oriented format

If a file does not begin with `{`, `InputParser` falls back to the legacy
line-oriented `key: value` parser. New configs should use strict JSON.

## Related documents

| Document | Content |
|----------|---------|
| [PARAMETERS.md](PARAMETERS.md) | All keys, defaults, units |
| [MECHANISMS.md](MECHANISMS.md) | Biological basis per Fix |
| [BATCH_RUNNER.md](BATCH_RUNNER.md) | Multi-run sweep manifests |
