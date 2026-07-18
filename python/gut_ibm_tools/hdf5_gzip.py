"""Whole-file gzip helpers for GutIBM HDF5 outputs."""

from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path


def gzip_hdf5_enabled_from_env(default: bool = False) -> bool:
    """Interpret GUTIBM_GZIP_HDF5 (true/false/1/0/yes/no/on/off)."""
    raw = os.environ.get("GUTIBM_GZIP_HDF5")
    if raw is None:
        return default
    value = raw.strip().lower()
    if value in {"1", "true", "yes", "on"}:
        return True
    if value in {"0", "false", "no", "off"}:
        return False
    return default


def gzip_hdf5_file(path: Path | str, *, compresslevel: int = 6) -> Path | None:
    """Replace an HDF5 file with a whole-file ``.gz`` sibling.

    Returns the ``.h5.gz`` path, or ``None`` when the input is missing / already
    gzipped. This is independent of HDF5-internal grid ``compression: gzip``.
    """
    hdf5_path = Path(path)
    if not hdf5_path.is_file():
        return None
    name = hdf5_path.name
    if name.endswith(".h5.gz") or hdf5_path.suffix == ".gz":
        return hdf5_path

    gz_path = Path(str(hdf5_path) + ".gz")
    before = hdf5_path.stat().st_size
    with hdf5_path.open("rb") as src, gzip.open(
        gz_path, "wb", compresslevel=compresslevel
    ) as dst:
        while True:
            chunk = src.read(1024 * 1024)
            if not chunk:
                break
            dst.write(chunk)
    hdf5_path.unlink()
    after = gz_path.stat().st_size
    print(
        f"gzipped HDF5 {hdf5_path.name}: {before} -> {after} bytes ({gz_path.name})",
        file=sys.stderr,
    )
    return gz_path


def maybe_gzip_hdf5_file(path: Path | str) -> Path | None:
    """Gzip ``path`` when ``GUTIBM_GZIP_HDF5`` is enabled."""
    if not gzip_hdf5_enabled_from_env(default=False):
        return None
    return gzip_hdf5_file(path)
