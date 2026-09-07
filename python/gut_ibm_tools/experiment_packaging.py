"""Shared fail-closed scaffolding for experiments/* deployment packages.

Deployment packages under ``experiments/`` generate and preflight JSON inputs
that must be attributable to an exact execution source commit and an immutable
container image.  These helpers centralize the repository-git wrapper, the
atomic manifest/input writer, the SHA-256 digest, and the deployment identity
checks so each package does not carry its own copy.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np

SHA40 = re.compile(r"[0-9a-f]{40}")
IMAGE_DIGEST = re.compile(r"[^\s@]+@sha256:[0-9a-f]{64}")


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *args], text=True, stderr=subprocess.STDOUT
    ).strip()


def write_json_atomic(path: Path, obj: dict) -> None:
    """Write ``obj`` as sorted, indented JSON atomically with an fsync."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(obj, indent=2, sort_keys=True) + "\n").encode("utf-8")
    handle, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_safe(value):
    """Convert a value for ``json.dumps(allow_nan=False)`` serialization."""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    return value


def require_immutable_identity(execution_sha: str, image_digest: str) -> None:
    """Refuse unless both identities are immutable deployment identities."""
    if not SHA40.fullmatch(execution_sha or ""):
        raise SystemExit(
            "REFUSED: --deployment requires --execution-source-sha as a full "
            "lowercase 40-hex commit"
        )
    if not IMAGE_DIGEST.fullmatch(image_digest or ""):
        raise SystemExit(
            "REFUSED: --deployment requires --image-digest <repo>@sha256:<64 hex>"
        )


def require_checkout_head(repo: Path, execution_sha: str) -> str:
    """Return HEAD, refusing unless it equals ``execution_sha`` in a git tree."""
    try:
        inside = git(repo, "rev-parse", "--is-inside-work-tree")
        head = git(repo, "rev-parse", "HEAD")
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(
            f"REFUSED: deployment generation requires a git checkout: {exc}"
        ) from exc
    if inside != "true" or head != execution_sha:
        raise SystemExit(
            f"REFUSED: --execution-source-sha {execution_sha} does not match "
            f"git HEAD {head}"
        )
    return head


def require_ancestor(repo: Path, ancestor: str, head: str) -> None:
    """Refuse unless ``ancestor`` is a commit ancestor of ``head``."""
    try:
        git(repo, "cat-file", "-e", ancestor + "^{commit}")
        ancestry = subprocess.call(
            ["git", "-C", str(repo), "merge-base", "--is-ancestor", ancestor, head],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(
            f"REFUSED: cannot verify required ancestry of {ancestor}: {exc}"
        ) from exc
    if ancestry != 0:
        raise SystemExit(
            f"REFUSED: required ancestor {ancestor} is not an ancestor of HEAD {head}"
        )


def require_package_clean(repo: Path, package_root: Path, package_files) -> None:
    """Refuse if any tracked package file differs from HEAD."""
    relative = package_root.relative_to(repo)
    dirty = git(
        repo,
        "status",
        "--porcelain",
        "--",
        *[str(relative / name) for name in package_files],
    )
    if dirty:
        raise SystemExit(
            "REFUSED: refinement package files differ from HEAD; commit them "
            "before deployment generation:\n" + dirty
        )
