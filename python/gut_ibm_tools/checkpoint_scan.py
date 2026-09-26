"""Discover and order GutIBM closed-restart checkpoint HDF5 files."""

from __future__ import annotations

import re
import warnings
from pathlib import Path

from .path_utils import PathValidationError, validate_path_syntax

STEP_H5_RE = re.compile(r"^step_(\d+)\.h5$", re.IGNORECASE)
CHECKPOINT_GLOB = "step_*.h5"


def validate_checkpoint_dir(path: str | Path) -> Path:
    """Validate that *path* is an existing directory (no ``..`` traversal)."""
    candidate = validate_path_syntax(path)
    if not candidate.exists():
        raise PathValidationError(f"checkpoint directory not found: {candidate}")
    if not candidate.is_dir():
        raise PathValidationError(f"checkpoint path is not a directory: {candidate}")
    return candidate.resolve()


def discover_checkpoints(directory: str | Path) -> list[Path]:
    """Return ``step_NNNNNN.h5`` paths sorted by embedded step index.

    Ignores incomplete upload artifacts such as ``step_000030.h5.6cF6bCBE``.
    """
    root = validate_checkpoint_dir(directory)
    found: list[tuple[int, Path]] = []
    for path in root.iterdir():
        if not path.is_file():
            continue
        match = STEP_H5_RE.fullmatch(path.name)
        if match is None:
            continue
        found.append((int(match.group(1)), path.resolve()))
    found.sort(key=lambda item: item[0])
    return [path for _, path in found]


def step_index_from_path(path: Path) -> int | None:
    """Parse the numeric step index from a checkpoint filename, if present."""
    match = STEP_H5_RE.fullmatch(path.name)
    if match is None:
        return None
    return int(match.group(1))


def warn_skip(path: Path, reason: str) -> None:
    """Emit a non-fatal warning for a skipped checkpoint."""
    warnings.warn(f"skipping checkpoint {path.name}: {reason}", stacklevel=2)
