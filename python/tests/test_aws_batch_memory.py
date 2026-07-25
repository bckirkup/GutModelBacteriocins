"""Golden + sensitivity tests for Batch memory-guard helpers."""

from __future__ import annotations

import pytest

from gut_ibm_tools.aws_batch_memory import (
    cgroup_free_mb,
    parse_mem_available_kb,
    sample_from_parts,
    should_stop_for_memory,
)

# Golden /proc/meminfo fragment (MemAvailable = 16777216 kB = 16384 MiB).
GOLDEN_MEMINFO = """\
MemTotal:       32768000 kB
MemFree:         8192000 kB
MemAvailable:   16777216 kB
Buffers:          100000 kB
"""


def test_parse_mem_available_golden() -> None:
    kb = parse_mem_available_kb(GOLDEN_MEMINFO)
    assert kb == 16_777_216
    sample = sample_from_parts(meminfo_text=GOLDEN_MEMINFO)
    assert sample.mem_available_mb == 16_384


def test_cgroup_free_golden_and_unlimited() -> None:
    # 28 GiB limit, 4 GiB used → 24 GiB free.
    free = cgroup_free_mb(28 * 1024**3, 4 * 1024**3)
    assert free == 24 * 1024
    assert cgroup_free_mb(2**63 - 1, 1024) is None


def test_should_stop_sensitivity_to_ram_threshold() -> None:
    sample = sample_from_parts(
        meminfo_text="MemAvailable: 2097152 kB\n",  # 2048 MiB
        gpu_free_mb=8000,
    )
    assert sample.mem_effective_free_mb == 2048
    # Boundary: free == min → do not stop; free < min → stop.
    assert should_stop_for_memory(sample, min_available_mb=2048) is False
    assert should_stop_for_memory(sample, min_available_mb=2049) is True


def test_should_stop_sensitivity_to_gpu_threshold() -> None:
    sample = sample_from_parts(
        meminfo_text=GOLDEN_MEMINFO,
        gpu_free_mb=400,
    )
    assert should_stop_for_memory(sample, min_gpu_free_mb=512) is True
    assert should_stop_for_memory(sample, min_gpu_free_mb=400) is False
    assert should_stop_for_memory(sample, min_gpu_free_mb=0) is False


def test_guard_disable_sensitivity() -> None:
    sample = sample_from_parts(meminfo_text="MemAvailable: 100 kB\n", gpu_free_mb=1)
    assert should_stop_for_memory(sample, guard_enabled=True) is True
    assert should_stop_for_memory(sample, guard_enabled=False) is False


def test_effective_prefers_tighter_cgroup() -> None:
    sample = sample_from_parts(
        meminfo_text=GOLDEN_MEMINFO,
        cgroup_max_bytes=8 * 1024**3,
        cgroup_current_bytes=6 * 1024**3,
    )
    assert sample.mem_available_mb == 16_384
    assert sample.mem_cgroup_free_mb == 2 * 1024
    assert sample.mem_effective_free_mb == pytest.approx(2 * 1024)
