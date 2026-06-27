"""Tests for gut_ibm_tools.path_utils."""

from __future__ import annotations

from pathlib import Path

import pytest

from gut_ibm_tools.path_utils import (
    PathValidationError,
    validate_input_path,
    validate_output_path,
    validate_path_syntax,
)


def test_validate_path_syntax_rejects_parent_traversal() -> None:
    with pytest.raises(PathValidationError, match="traversal"):
        validate_path_syntax("../etc/passwd")


def test_validate_path_syntax_rejects_null_byte() -> None:
    with pytest.raises(PathValidationError, match="null byte"):
        validate_path_syntax("safe\0path")


def test_validate_input_path_requires_existing_file(tmp_path: Path) -> None:
    target = tmp_path / "input.h5"
    target.write_bytes(b"data")
    resolved = validate_input_path(target)
    assert resolved == target.resolve()


def test_validate_output_path_rejects_symlink_in_world_writable_parent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    world_dir = tmp_path / "world_writable"
    world_dir.mkdir(mode=0o777)
    real_file = world_dir / "real.h5"
    real_file.write_bytes(b"data")
    link = world_dir / "linked.h5"
    link.symlink_to(real_file)

    monkeypatch.setattr(
        "gut_ibm_tools.path_utils._is_world_writable_directory",
        lambda path: path == world_dir.resolve(),
    )

    with pytest.raises(PathValidationError, match="symlink"):
        validate_output_path(link)


def test_validate_output_path_allows_private_temp_dir(tmp_path: Path) -> None:
    out = tmp_path / "nested" / "output.h5"
    out.parent.mkdir(parents=True)
    assert validate_output_path(out) == out
