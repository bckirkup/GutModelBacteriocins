---
name: sonarqube-python
description: SonarQube Python coding patterns for GutIBM tools. Use when writing or fixing python/gut_ibm_tools/ or python/tests/ to avoid float equality, path traversal, RNG, and complexity findings.
---

# SonarQube Python Patterns — GutIBM

Rule-keyed reference for `python/gut_ibm_tools/` and `python/tests/`.
Project: `bckirkup_GutModelBacteriocins`.

See `.agents/skills/sonarqube-gutibm/SKILL.md` for batch PR workflow.

## Float equality (S1244) — BUG, not smell

Never use `==` or `!=` on float literals or computed floats in tests or assertions:

```python
# Bad — Sonar BUG
assert metrics["score"] == 0.0
assert result["fraction"] == 1.0

# Good — exact literals
assert metrics["score"] == pytest.approx(0.0)
assert result["fraction"] == pytest.approx(1.0)

# Good — computed / arrays
assert value == pytest.approx(2 / 3, rel=1e-3)
np.testing.assert_allclose(times, [0.0, 3600.0])
```

Existing patterns: `test_validation.py`, `test_hdf5_writer_roundtrip.py`, `test_eari_vadi_regression.py`.

## Path safety (pythonsecurity:S8707) — VULNERABILITY

User-supplied paths (CLI args, env) must be validated **before** any filesystem call (`mkdir`, `open`, `Path.read_text`).

Use `gut_ibm_tools.path_utils`:

| Function | Use |
|----------|-----|
| `validate_path_syntax()` | Reject `..`, null bytes, empty paths |
| `validate_input_path()` | Read existing files; returns resolved path |
| `validate_output_path()` | Write when parent already exists |
| `prepare_output_file()` | Validate then `mkdir -p` then re-validate — **use for new output files** |

```python
from gut_ibm_tools.path_utils import prepare_output_file, validate_input_path

def write_golden(metrics: dict, path: str | Path) -> None:
    out = prepare_output_file(path)          # validates BEFORE mkdir
    with open(out, "w", encoding="utf-8") as f:
        json.dump(metrics, f)

def load_golden(path: str | Path) -> dict:
    safe_path = validate_input_path(path)
    with open(safe_path, encoding="utf-8") as f:
        return json.load(f)
```

**Anti-pattern** (triggers S8707):

```python
out = validate_path_syntax(path)
out.parent.mkdir(parents=True, exist_ok=True)  # sink before full validation
validate_output_path(out)
```

Reference: `python/gut_ibm_tools/validation_regression.py` (`write_golden`, `write_fish_golden`).

Tests: `python/tests/test_path_utils.py`.

## Random number generation (S6709, S6711)

Use NumPy Generator API — not legacy global RNG:

```python
rng = np.random.default_rng(seed)
indices = rng.choice(n, m, replace=False)
values = rng.uniform(lo, hi, size=(m, 3))
```

Pass `rng` as an optional parameter; default to `np.random.default_rng()` when `None`.

## String literal duplication (S1192) — CRITICAL severity smell

Extract repeated literals (issue references, metric names) to module constants:

```python
ISSUE_25_REF = "issue #25"

FISH_TARGETS = {
    "immunity_hcr": {"references": [ISSUE_25_REF, "VADI"]},
}
```

Hotspot: `python/gut_ibm_tools/validation_regression.py`.

## Cognitive complexity (S3776) — CRITICAL severity smell

Keep functions ≤ 15 cognitive complexity. For CLI `main()`:

- Parse args in `main()`
- Delegate to `run_validation()`, `_write_golden_outputs()`, `_print_failures()`
- One `if` branch per subcommand handler

## Naming (S117)

Follow PEP 8: `snake_case` for functions and variables.

## Lint and test gates

Every Python change:

```bash
cd python && ruff check .
cd python && pytest tests/ -v -m "not integration"
```

CI job `python-lint` runs the same checks on every PR.

## New Python file checklist

```
□ Float assertions use pytest.approx or np.isclose / assert_allclose
□ User paths go through path_utils (prepare_output_file for writes)
□ Randomness uses np.random.default_rng()
□ Repeated string literals extracted to constants
□ CLI logic split into helpers if complexity > 15
□ ruff check passes
```

## Reference files

- Path utils: `python/gut_ibm_tools/path_utils.py`
- Golden writers: `python/gut_ibm_tools/validation_regression.py`
- Path tests: `python/tests/test_path_utils.py`
- Float assertion examples: `python/tests/test_hdf5_writer_roundtrip.py`
