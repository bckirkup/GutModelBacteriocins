"""Filesystem path validation helpers for GutIBM Python tools."""

from __future__ import annotations

import os
import stat
from pathlib import Path


class PathValidationError(ValueError):
    """Raised when a constructed path fails safety checks."""


def _path_has_parent_traversal(path: Path) -> bool:
    return ".." in path.parts


def validate_path_syntax(path: str | Path) -> Path:
    """Reject empty paths, null bytes, and parent-directory traversal."""
    if isinstance(path, str) and "\0" in path:
        raise PathValidationError("path contains null byte")

    resolved = Path(path)
    if not resolved.parts:
        raise PathValidationError("empty path")
    if _path_has_parent_traversal(resolved):
        raise PathValidationError("path contains parent-directory traversal ('..')")
    return resolved


def _is_world_writable_directory(path: Path) -> bool:
    try:
        mode = path.stat(follow_symlinks=False).st_mode
    except OSError:
        return False
    return stat.S_ISDIR(mode) and bool(mode & stat.S_IWOTH)


def _reject_symlink_in_world_writable_parent(path: Path, operation: str) -> None:
    parent = path.parent if path.parent.parts else Path(".")
    if not parent.exists():
        return

    parent_resolved = parent.resolve()
    if not _is_world_writable_directory(parent_resolved):
        return

    if parent.is_symlink():
        raise PathValidationError(
            f"refusing to {operation} via symlinked parent in world-writable directory: {parent}"
        )

    if path.exists() and path.is_symlink():
        raise PathValidationError(
            f"refusing to {operation} through symlink in world-writable directory: {path}"
        )


def validate_input_path(path: str | Path) -> Path:
    """Validate a readable file path and return its canonical form."""
    candidate = validate_path_syntax(path)
    if not candidate.exists():
        raise PathValidationError(f"input file not found: {candidate}")

    _reject_symlink_in_world_writable_parent(candidate, "read")

    if not candidate.is_file():
        raise PathValidationError(f"input path is not a regular file: {candidate}")

    return candidate.resolve()


def validate_output_path(path: str | Path) -> Path:
    """Validate a writable file path before creating or overwriting it."""
    candidate = validate_path_syntax(path)
    parent = candidate.parent if candidate.parent.parts else Path(".")

    if not parent.exists():
        raise PathValidationError(f"output directory does not exist: {parent}")
    if not parent.is_dir():
        raise PathValidationError(f"output parent is not a directory: {parent}")

    _reject_symlink_in_world_writable_parent(candidate, "write")

    if candidate.exists():
        if candidate.is_symlink():
            raise PathValidationError(f"refusing to overwrite symlink: {candidate}")
        if not candidate.is_file():
            raise PathValidationError(
                f"output path exists and is not a regular file: {candidate}"
            )

    return candidate
