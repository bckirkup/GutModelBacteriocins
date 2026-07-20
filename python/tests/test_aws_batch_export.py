"""Tests for the AWS Batch array-export helper (Phase 3).

Parity anchor: expanded index count and per-index input contents must match
``batch_runner`` expansion (same ``parse_batch_config`` / ``build_job_config``).
Sensitivity probe: sweep values must actually differ across the array (guards
against dead-wiring where every input is identical). S3 / submit are mocked.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gut_ibm_tools.aws_batch_export import (
    ENV_CHECKPOINT_PREFIX,
    ENV_INPUT_PREFIX,
    ENV_MPI_RANKS,
    ENV_OUTPUT_PREFIX,
    export_array,
)
from gut_ibm_tools.batch_config import (
    BatchConfigError,
    apply_overrides,
    load_simulation_template,
    parse_batch_config,
    strip_metadata_keys,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_DIR = "experiments/diversity_campaign/stage3_campaign"
KD_SWEEP = f"{CAMPAIGN_DIR}/batch_kd_sweep.json"
BASELINE = f"{CAMPAIGN_DIR}/batch_baseline.json"

INPUT_PREFIX = "s3://bucket/campaign/jobs"
OUTPUT_PREFIX = "s3://bucket/campaign/out"
CHECKPOINT_PREFIX = "s3://bucket/campaign/ckpt"
QUEUE = "gutibm-gpu-campaign"
JOB_DEF = "gutibm-cuda-campaign"


class _Recorder:
    def __init__(self, submit_stdout: str = '{"jobId": "abc-123"}') -> None:
        self.uploads: list[tuple[str, str]] = []
        self.submits: list[list[str]] = []
        self._submit_stdout = submit_stdout

    def s3_put(self, content: str, s3_uri: str) -> None:
        self.uploads.append((content, s3_uri))

    def submit(self, argv: list[str]) -> str:
        self.submits.append(argv)
        return self._submit_stdout


def _run(batch: str, recorder: _Recorder, *, dry_run: bool = False, checkpoint: str | None = CHECKPOINT_PREFIX):
    return export_array(
        REPO_ROOT / batch,
        input_prefix=INPUT_PREFIX,
        output_prefix=OUTPUT_PREFIX,
        checkpoint_prefix=checkpoint,
        job_queue=QUEUE,
        job_definition=JOB_DEF,
        dry_run=dry_run,
        repo_root=REPO_ROOT,
        s3_put=recorder.s3_put,
        submit=recorder.submit,
    )


def _runner_configs(batch: str) -> list[dict]:
    """Rebuild per-index configs the way batch_runner would (independent path)."""
    settings, jobs = parse_batch_config(
        REPO_ROOT / batch, repo_root=REPO_ROOT, require_binary=False,
    )
    template = strip_metadata_keys(load_simulation_template(settings.base_config))
    configs = []
    for job in jobs:
        cfg = json.loads(json.dumps(template))
        apply_overrides(cfg, job.overrides)
        configs.append(cfg)
    return configs


@pytest.mark.parametrize(("batch", "expected_size"), [(KD_SWEEP, 12), (BASELINE, 3)])
def test_array_size_matches_batch_runner(batch: str, expected_size: int) -> None:
    recorder = _Recorder()
    result = _run(batch, recorder, dry_run=True)
    _settings, jobs = parse_batch_config(
        REPO_ROOT / batch, repo_root=REPO_ROOT, require_binary=False,
    )
    assert result.array_size == expected_size
    assert len(jobs) == expected_size
    assert len(result.configs) == expected_size


@pytest.mark.parametrize("batch", [KD_SWEEP, BASELINE])
def test_per_index_config_matches_batch_runner(batch: str) -> None:
    recorder = _Recorder()
    result = _run(batch, recorder, dry_run=True)
    expected = _runner_configs(batch)
    assert len(result.configs) == len(expected)
    for exported, ref in zip(result.configs, expected, strict=True):
        # Overrides propagate identically to batch_runner; export adds gpu_enabled
        # and a stable hdf5_file placeholder.
        ref_with_gpu = {**ref, "gpu_enabled": True, "hdf5_file": "output.h5"}
        assert exported["gpu_enabled"] is True
        assert exported["hdf5_file"] == "output.h5"
        assert json.dumps(exported, sort_keys=True) == json.dumps(
            ref_with_gpu, sort_keys=True,
        )


def test_kd_sweep_golden_expansion_order() -> None:
    recorder = _Recorder()
    result = _run(KD_SWEEP, recorder, dry_run=True)
    # sorted sweep keys => kd_corrinoid_btuB (outer) x seed (inner).
    assert result.configs[0]["seed"] == 42
    assert result.configs[0]["kd_corrinoid_btuB"] == pytest.approx(1e-9)
    assert result.configs[1]["seed"] == 43
    assert result.configs[-1]["seed"] == 44
    assert result.configs[-1]["kd_corrinoid_btuB"] == pytest.approx(1e-6)


def test_sweep_values_differ_across_array() -> None:
    """Sensitivity: sweep must actually vary inputs, not emit identical configs."""
    recorder = _Recorder()
    result = _run(KD_SWEEP, recorder, dry_run=True)
    kds = {cfg["kd_corrinoid_btuB"] for cfg in result.configs}
    seeds = {cfg["seed"] for cfg in result.configs}
    assert len(kds) == 4
    assert seeds == {42, 43, 44}
    fingerprints = {json.dumps(cfg, sort_keys=True) for cfg in result.configs}
    assert len(fingerprints) == len(result.configs)


def test_upload_uris_and_submit_overrides() -> None:
    recorder = _Recorder()
    result = _run(KD_SWEEP, recorder)
    assert [uri for _content, uri in recorder.uploads] == [
        f"{INPUT_PREFIX}/{i}/input.json" for i in range(12)
    ]
    assert result.input_uris == [f"{INPUT_PREFIX}/{i}/input.json" for i in range(12)]
    assert len(recorder.submits) == 1
    argv = recorder.submits[0]
    assert "--array-properties" in argv
    assert argv[argv.index("--array-properties") + 1] == "size=12"
    overrides = json.loads(argv[argv.index("--container-overrides") + 1])
    env = {item["name"]: item["value"] for item in overrides["environment"]}
    assert env[ENV_INPUT_PREFIX] == INPUT_PREFIX
    assert env[ENV_OUTPUT_PREFIX] == OUTPUT_PREFIX
    assert env[ENV_CHECKPOINT_PREFIX] == CHECKPOINT_PREFIX
    assert env[ENV_MPI_RANKS] == "1"
    assert result.job_id == "abc-123"


def test_uploaded_content_round_trips() -> None:
    recorder = _Recorder()
    result = _run(BASELINE, recorder)
    for (content, _uri), cfg in zip(recorder.uploads, result.configs, strict=True):
        assert json.loads(content) == cfg


def test_dry_run_uploads_nothing() -> None:
    recorder = _Recorder()
    result = _run(KD_SWEEP, recorder, dry_run=True)
    assert recorder.uploads == []
    assert recorder.submits == []
    assert result.job_id is None


def test_checkpoint_prefix_optional() -> None:
    recorder = _Recorder()
    _run(BASELINE, recorder, checkpoint=None)
    argv = recorder.submits[0]
    overrides = json.loads(argv[argv.index("--container-overrides") + 1])
    names = {item["name"] for item in overrides["environment"]}
    assert ENV_CHECKPOINT_PREFIX not in names


def test_non_s3_prefix_rejected() -> None:
    recorder = _Recorder()
    with pytest.raises(BatchConfigError):
        export_array(
            REPO_ROOT / BASELINE,
            input_prefix="/local/path",
            output_prefix=OUTPUT_PREFIX,
            job_queue=QUEUE,
            job_definition=JOB_DEF,
            repo_root=REPO_ROOT,
            s3_put=recorder.s3_put,
            submit=recorder.submit,
        )
