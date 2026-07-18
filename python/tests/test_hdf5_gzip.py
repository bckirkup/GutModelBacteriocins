"""Tests for whole-file HDF5 gzip helpers."""

from __future__ import annotations

import gzip
from pathlib import Path

from gut_ibm_tools.hdf5_gzip import (
    gzip_hdf5_enabled_from_env,
    gzip_hdf5_file,
    maybe_gzip_hdf5_file,
)


def test_gzip_hdf5_file_replaces_with_gz(tmp_path: Path, monkeypatch) -> None:
    hdf5 = tmp_path / "output.h5"
    payload = b"HDF5-fake-content-" + (b"\0" * 4096)
    hdf5.write_bytes(payload)

    gz_path = gzip_hdf5_file(hdf5)
    assert gz_path == tmp_path / "output.h5.gz"
    assert gz_path is not None
    assert gz_path.is_file()
    assert not hdf5.exists()
    assert gz_path.stat().st_size < len(payload)
    assert gzip.decompress(gz_path.read_bytes()) == payload


def test_gzip_hdf5_file_noop_when_missing(tmp_path: Path) -> None:
    assert gzip_hdf5_file(tmp_path / "missing.h5") is None


def test_maybe_gzip_respects_env(tmp_path: Path, monkeypatch) -> None:
    hdf5 = tmp_path / "run.h5"
    hdf5.write_bytes(b"abc" * 1000)
    monkeypatch.delenv("GUTIBM_GZIP_HDF5", raising=False)
    assert maybe_gzip_hdf5_file(hdf5) is None
    assert hdf5.is_file()

    monkeypatch.setenv("GUTIBM_GZIP_HDF5", "true")
    assert gzip_hdf5_enabled_from_env() is True
    gz = maybe_gzip_hdf5_file(hdf5)
    assert gz is not None
    assert gz.name == "run.h5.gz"
    assert not hdf5.exists()
