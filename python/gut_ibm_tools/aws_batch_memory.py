"""Pure helpers for Batch memory-guard logic (mirrors deploy/aws/entry.sh).

Kept in Python so CI can golden-test MemAvailable parsing and threshold
sensitivity without running the container entrypoint.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MemorySample:
    mem_available_mb: int | None
    mem_cgroup_free_mb: int | None
    mem_effective_free_mb: int | None
    gpu_free_mb: int | None


def parse_mem_available_kb(meminfo_text: str) -> int | None:
    """Parse MemAvailable kilobytes from a ``/proc/meminfo`` dump."""
    for line in meminfo_text.splitlines():
        if line.startswith("MemAvailable:"):
            parts = line.split()
            if len(parts) >= 2 and parts[1].isdigit():
                return int(parts[1])
            return None
    return None


def cgroup_free_mb(max_bytes: int | None, current_bytes: int | None) -> int | None:
    """Free MiB under a cgroup memory.max / memory.current pair."""
    if max_bytes is None or current_bytes is None:
        return None
    # Treat huge "unlimited" legacy caps as unset (same idea as entry.sh).
    if max_bytes >= 1_000_000_000_000:
        return None
    free = max_bytes - current_bytes
    if free < 0:
        free = 0
    return free // (1024 * 1024)


def effective_free_mb(
    host_available_mb: int | None,
    cgroup_free: int | None,
) -> int | None:
    """Min of host MemAvailable and cgroup free when both are known."""
    if host_available_mb is not None and cgroup_free is not None:
        return min(host_available_mb, cgroup_free)
    if cgroup_free is not None:
        return cgroup_free
    return host_available_mb


def sample_from_parts(
    *,
    meminfo_text: str | None = None,
    cgroup_max_bytes: int | None = None,
    cgroup_current_bytes: int | None = None,
    gpu_free_mb: int | None = None,
) -> MemorySample:
    """Build a MemorySample from injectable strings/ints (unit-test friendly)."""
    host_mb: int | None = None
    if meminfo_text is not None:
        kb = parse_mem_available_kb(meminfo_text)
        host_mb = None if kb is None else kb // 1024
    cg = cgroup_free_mb(cgroup_max_bytes, cgroup_current_bytes)
    return MemorySample(
        mem_available_mb=host_mb,
        mem_cgroup_free_mb=cg,
        mem_effective_free_mb=effective_free_mb(host_mb, cg),
        gpu_free_mb=gpu_free_mb,
    )


def should_stop_for_memory(
    sample: MemorySample,
    *,
    min_available_mb: int = 2048,
    min_gpu_free_mb: int = 512,
    guard_enabled: bool = True,
) -> bool:
    """True when the entrypoint memory guard should SIGTERM gut_ibm."""
    if not guard_enabled:
        return False
    eff = sample.mem_effective_free_mb
    if eff is not None and eff < min_available_mb:
        return True
    if min_gpu_free_mb > 0 and sample.gpu_free_mb is not None:
        if sample.gpu_free_mb < min_gpu_free_mb:
            return True
    return False
