"""Tests for AWS Batch array output QA (summary fingerprint gate)."""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest
from gut_ibm_tools.aws_batch_qa import (
    format_report,
    list_array_indices,
    qa_array,
    summarize_hdf5,
)


def _write_summary_h5(
    path: Path,
    *,
    n_agents: int,
    boundary_deaths: int,
    washout_deaths: int = 0,
    mean_carbon: float = 0.002,
    with_agents: bool = False,
) -> None:
    with h5py.File(path, "w") as f:
        f.attrs["gutibm_version"] = 4
        f.attrs["nx"] = 20
        f.attrs["ny"] = 20
        f.attrs["nz"] = 10
        summary = f.require_group("summary").require_group("step_000002")
        summary.create_dataset("time", data=np.array(120.0))
        summary.create_dataset("step", data=np.array(2, dtype=np.int32))
        summary.create_dataset("n_total", data=np.array(n_agents, dtype=np.int32))
        summary.create_dataset("num_agents", data=np.array(n_agents, dtype=np.int32))
        events = summary.require_group("events")
        events.create_dataset("boundary_deaths", data=np.array(boundary_deaths, dtype=np.int32))
        events.create_dataset("washout_deaths", data=np.array(washout_deaths, dtype=np.int32))
        events.create_dataset("colicin_kills", data=np.array(0, dtype=np.int32))
        chem = summary.require_group("chem")
        chem.create_dataset("mean_carbon", data=np.array(mean_carbon))
        if with_agents:
            agents = f.require_group("agents").require_group("step_000002")
            agents.create_dataset("id", data=np.arange(n_agents, dtype=np.int64))
            agents.create_dataset("x", data=np.zeros(n_agents))
            agents.create_dataset("y", data=np.zeros(n_agents))
            agents.create_dataset("z", data=np.zeros(n_agents))


def _gzip_file(src: Path, dest: Path) -> None:
    with src.open("rb") as fh_in, gzip.open(dest, "wb") as fh_out:
        shutil.copyfileobj(fh_in, fh_out)


class _FakeStore:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.downloads: list[tuple[str, Path]] = []

    def aws_run(self, argv: list[str]) -> str:
        assert argv[:3] == ["aws", "s3", "ls"]
        assert argv[3].startswith("s3://")
        lines = []
        for child in sorted(self.root.iterdir()):
            if child.is_dir() and child.name.isdigit():
                lines.append(f"PRE {child.name}/")
        return "\n".join(lines) + "\n"

    def s3_download(self, s3_uri: str, dest: Path) -> None:
        self.downloads.append((s3_uri, dest))
        # s3://bucket/out/0/output.h5.gz -> root/0/output.h5.gz
        parts = s3_uri.split("/")
        index = parts[-2]
        name = parts[-1]
        src = self.root / index / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)


def test_summarize_hdf5_fingerprint_changes_with_agents(tmp_path: Path) -> None:
    a = tmp_path / "a.h5"
    b = tmp_path / "b.h5"
    _write_summary_h5(a, n_agents=13, boundary_deaths=7)
    _write_summary_h5(b, n_agents=12, boundary_deaths=7)
    sa = summarize_hdf5(a, index=0, seed=4092)
    sb = summarize_hdf5(b, index=1, seed=4093)
    assert sa.final_agents == 13
    assert sb.final_agents == 12
    assert sa.fingerprint != sb.fingerprint
    assert sa.has_agents_layer is False


def test_qa_array_distinct_seeds_and_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    store_root = tmp_path / "s3"
    for index, seed, agents, deaths in (
        (0, 4092, 13, 7),
        (1, 4093, 12, 7),
    ):
        child = store_root / str(index)
        child.mkdir(parents=True)
        h5 = child / "raw.h5"
        _write_summary_h5(h5, n_agents=agents, boundary_deaths=deaths)
        _gzip_file(h5, child / "output.h5.gz")
        (child / "input.json").write_text(
            json.dumps({"seed": seed, "gpu_enabled": True}),
            encoding="utf-8",
        )

    store = _FakeStore(store_root)
    work = Path("work")
    report = qa_array(
        output_prefix="s3://bucket/out",
        input_prefix="s3://bucket/jobs",
        work_dir=work,
        aws_run=store.aws_run,
        s3_download=store.s3_download,
    )
    assert report.distinct_fingerprints is True
    assert report.distinct_seeds is True
    assert [row.seed for row in report.rows] == [4092, 4093]
    assert [row.final_agents for row in report.rows] == [13, 12]
    text = format_report(report)
    assert "4092" in text
    assert "4093" in text
    assert "distinct_fingerprints=True" in text
    assert any("summary-only" in w for w in report.warnings)


def test_qa_array_fails_when_identical(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    store_root = tmp_path / "s3"
    for index in (0, 1):
        child = store_root / str(index)
        child.mkdir(parents=True)
        h5 = child / "raw.h5"
        _write_summary_h5(h5, n_agents=13, boundary_deaths=7)
        _gzip_file(h5, child / "output.h5.gz")
        (child / "input.json").write_text(json.dumps({"seed": 4092}), encoding="utf-8")

    store = _FakeStore(store_root)
    with pytest.raises(AssertionError, match="not distinct"):
        qa_array(
            output_prefix="s3://bucket/out",
            input_prefix="s3://bucket/jobs",
            work_dir=Path("work"),
            aws_run=store.aws_run,
            s3_download=store.s3_download,
        )


def test_list_array_indices_parses_pre_lines() -> None:
    def aws_run(argv: list[str]) -> str:
        return "                           PRE 0/\n                           PRE 1/\n"

    assert list_array_indices("s3://bucket/out", aws_run=aws_run) == [0, 1]
