# GutIBM Input Config Format

GutIBM simulation configs are **strict JSON objects**. They validate with standard
JSON tools (`python -m json.tool`, `jq`, IDE parsers) and are loaded by
`InputParser::parse()` via the built-in JSON document parser.

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
| `total_time`, `bio_dt`, … | number | See [PARAMETERS.md](PARAMETERS.md) |
| `initial_strains` | array | Strain objects with `type`, `count`, `plasmids`, … |
| `fixes` | string[] | Fix plugin names in execution order |

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
  "hdf5_every": 60
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
