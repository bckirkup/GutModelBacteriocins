"""Filesystem path validation helpers for GutIBM Python tools."""

from __future__ import annotations

import json
import re
import stat
from pathlib import Path
from typing import Any


class PathValidationError(ValueError):
    """Raised when a constructed path fails safety checks."""


_SAFE_PATH_SEGMENT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_PARENT_TRAVERSAL_MSG = "path contains parent-directory traversal ('..')"


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
        raise PathValidationError(_PARENT_TRAVERSAL_MSG)
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


def _output_parent_directory(parent: Path) -> Path:
    return parent if parent.parts else Path(".")


def _validate_output_parent_directory(parent: Path) -> None:
    if parent.exists():
        if not parent.is_dir():
            raise PathValidationError(f"output parent is not a directory: {parent}")
        return

    ancestor = parent
    while not ancestor.exists() and ancestor.parts:
        ancestor = _output_parent_directory(ancestor.parent)
    if ancestor.exists() and not ancestor.is_dir():
        raise PathValidationError(f"output ancestor is not a directory: {ancestor}")


def _validate_existing_output_file(candidate: Path) -> None:
    if not candidate.exists():
        return
    if candidate.is_symlink():
        raise PathValidationError(f"refusing to overwrite symlink: {candidate}")
    if not candidate.is_file():
        raise PathValidationError(
            f"output path exists and is not a regular file: {candidate}"
        )


def _validate_output_parent(candidate: Path) -> None:
    parent = _output_parent_directory(candidate.parent)
    _validate_output_parent_directory(parent)
    _reject_symlink_in_world_writable_parent(candidate, "write")
    _validate_existing_output_file(candidate)


def validate_output_path(path: str | Path) -> Path:
    """Validate a writable file path before creating or overwriting it."""
    candidate = validate_path_syntax(path)
    parent = candidate.parent if candidate.parent.parts else Path(".")

    if not parent.exists():
        raise PathValidationError(f"output directory does not exist: {parent}")
    if not parent.is_dir():
        raise PathValidationError(f"output parent is not a directory: {parent}")

    _validate_output_parent(candidate)
    return candidate


def _ensure_output_within_cwd(candidate: Path) -> Path:
    """Reject output paths that resolve outside the current working directory."""
    cwd = Path.cwd().resolve()
    resolved = candidate.resolve()
    if not resolved.is_relative_to(cwd):
        raise PathValidationError(
            f"output path must resolve under the working directory: {candidate}"
        )
    return candidate


def _sanitize_path_segment(segment: str) -> str:
    """Allowlist a single path component to break S8707 taint from CLI arguments."""
    if segment in {".", ".."}:
        raise PathValidationError(_PARENT_TRAVERSAL_MSG)
    if not _SAFE_PATH_SEGMENT_RE.fullmatch(segment):
        raise PathValidationError(f"unsafe path segment: {segment!r}")
    return segment


def _trusted_output_path(candidate: Path) -> Path:
    """Rebuild an output path from cwd + sanitized segments (pythonsecurity:S8707)."""
    cwd = Path.cwd().resolve()
    resolved = candidate.resolve()
    if not resolved.is_relative_to(cwd):
        raise PathValidationError(
            f"output path must resolve under the working directory: {candidate}"
        )

    trusted = cwd
    for segment in resolved.relative_to(cwd).parts:
        trusted = trusted / _sanitize_path_segment(segment)

    if trusted.resolve() != resolved:
        raise PathValidationError(f"output path failed trust rebuild: {candidate}")
    return trusted


def _mkdir_validated_parents(parent: Path) -> None:
    """Create missing parent directories after the full path has been validated."""
    if parent.exists():
        if not parent.is_dir():
            raise PathValidationError(f"output parent is not a directory: {parent}")
        return

    to_create: list[Path] = []
    current = parent
    while not current.exists() and current.parts:
        if _path_has_parent_traversal(current):
            raise PathValidationError(_PARENT_TRAVERSAL_MSG)
        to_create.append(current)
        current = current.parent if current.parent.parts else Path(".")

    if current.exists() and not current.is_dir():
        raise PathValidationError(f"output ancestor is not a directory: {current}")

    for directory in reversed(to_create):
        directory.mkdir(exist_ok=True)


def prepare_output_file(path: str | Path) -> Path:
    """Validate an output file path and create parent directories if needed."""
    candidate = validate_path_syntax(path)
    parent = candidate.parent if candidate.parent.parts else Path(".")

    _validate_output_parent(candidate)
    _validate_output_parent_directory(parent)
    _reject_symlink_in_world_writable_parent(candidate, "write")
    _validate_existing_output_file(candidate)
    _mkdir_validated_parents(parent)
    return validate_output_path(candidate)


def write_text_file(path: str | Path, text: str) -> None:
    """Write text to a validated output path (must resolve under cwd)."""
    candidate = validate_path_syntax(path)
    candidate = _ensure_output_within_cwd(candidate)
    candidate = prepare_output_file(candidate)
    trusted_path = _trusted_output_path(candidate)
    _validate_output_parent(trusted_path)
    with trusted_path.open("w", encoding="utf-8") as handle:
        handle.write(text)


def write_json_file(path: str | Path, payload: Any, *, indent: int = 2) -> None:
    """Write JSON to a validated output path."""
    write_text_file(path, json.dumps(payload, indent=indent) + "\n")
