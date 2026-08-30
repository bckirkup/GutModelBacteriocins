#!/usr/bin/env python3
"""Bulk-resolve accepted SonarCloud technical-debt smells as Won't Fix.

SonarCloud automatic analysis ignores sonar.issue.ignore.multicriteria in
sonar-project.properties. After mechanical cleanup merges, run this script
with a project-admin token to clear the dashboard for accepted debt rules.

Usage:
  export SONAR_TOKEN=...   # SonarCloud user token with Administer Issues
  python3 scripts/sonar_wont_fix_debt.py [--dry-run]

See docs/SONARQUBE_PLAN.md.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

PROJECT = "bckirkup_GutModelBacteriocins"
BASE = "https://sonarcloud.io/api"

# Accepted debt, per rule, with the reason recorded on each resolution.
# Every other open rule is fixed in code — do not add a rule here to avoid
# doing the work, and never add a security rule.
DEBT_RULES: dict[str, str] = {
    # Language standard: the fix needs a C++23 *library* feature, gated on
    # __cplusplus > 202002L, so it does not compile at CMAKE_CXX_STANDARD 20
    # even under GCC 13 (the compiler the project actually builds with).
    # Bumping the standard is not a cleanup: CMAKE_CUDA_STANDARD is 17, and
    # shared __host__ __device__ headers are compiled by nvcc.
    "cpp:S7034": "std::string::contains is C++23; project is C++20",
    "cpp:S7035": "std::to_underlying is C++23; project is C++20",
    # Numerical reproducibility: std::lerp/std::midpoint are not bit-identical
    # to the current arithmetic. These sites are the Robin correction-table
    # interpolation and the metabolic-mode blend, whose outputs are compared
    # against Python oracles at ~1e-9 and are regression-guarded; changing the
    # arithmetic for a style rule trades reproducibility for nothing.
    "cpp:S6179": "std::lerp is not bit-identical here; FP reproducibility",
    # cpp:S6185 and cpp:S6484 (std::format) are NOT debt: <format> is a C++20
    # header, libstdc++ 13 has it, and src/io/hdf5_writer.cpp already uses it.
    # They are fixed in code. Do not re-add them on a GCC-11 assumption.
    # cpp:S8379 is deliberately NOT listed yet. Most of its findings are
    # already synchronized by a mechanism the rule cannot see (OpenMP atomic
    # updates and per-thread slots for the Green's-function diagnostics,
    # serial-only mutation for the Fix vector and the HDF5 provenance flag),
    # but triaging it found one genuine race, fixed in PR #371, and the two
    # FixMetabolism findings are still under audit. Resolving the family
    # wholesale would resolve those two as well. Add it here only once the
    # audit lands.
    # Complexity and architecture of a research prototype: NUFEB-style Fix
    # plugins, the diffusion kernels, and the config parser. Refactoring these
    # is a redesign, not a cleanup, and would put the scientific code at risk.
    "cpp:S107": "diffusion/GPU APIs need a context-struct redesign",
    "cpp:S134": "nesting in hot kernels, receptor, and GPU paths",
    "cpp:S3776": "parser/HDF5/GPU complexity; redesign, not cleanup",
    "python:S3776": "batch CLI and analysis complexity; redesign",
    "cpp:S1820": "Simulation/GPU type size is the architecture",
    "cpp:S1448": "Simulation/GPU method count is the architecture",
    "cpp:S995": "GPU buffer mutability",
    "cpp:S5008": "void* is the HDF5 C API buffer type",
    "cpp:S3656": "protected members are the NUFEB-style Fix base contract",
    "cpp:S924": "nested break is coupled to Simulation::run",
}

COMMENT_PREFIX = "Accepted GutIBM debt — see docs/SONARQUBE_PLAN.md. "

REOPEN_COMMENT = (
    "Reopened: the reason recorded on this resolution was wrong. "
    "See docs/SONARQUBE_PLAN.md — this rule is fixed in code."
)


def _request(
    method: str,
    path: str,
    token: str,
    params: dict[str, Any] | None = None,
    data: dict[str, Any] | None = None,
) -> dict[str, Any]:
    query = urllib.parse.urlencode(params or {}, doseq=True)
    url = f"{BASE}{path}"
    if query:
        url = f"{url}?{query}"
    body = None
    headers = {"Authorization": f"Bearer {token}"}
    if data is not None:
        body = urllib.parse.urlencode(data, doseq=True).encode()
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            raw = resp.read().decode()
            if not raw:
                return {}
            return json.loads(raw)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")
        raise SystemExit(f"HTTP {exc.code} on {path}: {detail}") from exc


def list_open_debt_issues(
    token: str, rules: tuple[str, ...] = ()
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    page = 1
    selected = rules or tuple(DEBT_RULES)
    while True:
        payload = _request(
            "GET",
            "/issues/search",
            token,
            params={
                "componentKeys": PROJECT,
                "resolved": "false",
                "rules": ",".join(selected),
                "ps": 100,
                "p": page,
            },
        )
        batch = payload.get("issues", [])
        issues.extend(batch)
        total = payload.get("total", 0)
        if page * 100 >= total or not batch:
            break
        page += 1
    return issues


def reopen_rules(token: str, rules: tuple[str, ...], dry_run: bool) -> int:
    keys: list[str] = []
    page = 1
    while True:
        payload = _request(
            "GET",
            "/issues/search",
            token,
            params={
                "componentKeys": PROJECT,
                "resolutions": "WONTFIX",
                "rules": ",".join(rules),
                "ps": 100,
                "p": page,
            },
        )
        batch = payload.get("issues", [])
        keys.extend(issue["key"] for issue in batch)
        if page * 100 >= int(payload.get("total", 0)) or not batch:
            break
        page += 1

    print(f"Won't Fix issues to reopen for {','.join(rules)}: {len(keys)}")
    for start in range(0, len(keys), 400):
        chunk = keys[start:start + 400]
        if dry_run:
            print(f"[dry-run] REOPEN {len(chunk)} issue(s)")
            continue
        _request(
            "POST",
            "/issues/bulk_change",
            token,
            data={
                "issues": ",".join(chunk),
                "do_transition": "reopen",
                "comment": REOPEN_COMMENT,
            },
        )
        print(f"REOPEN {len(chunk)} issue(s)")
    return 0


def bulk_wont_fix(
    token: str, issue_keys: list[str], comment: str, dry_run: bool
) -> None:
    # SonarCloud bulk_change accepts up to 500 keys per call.
    chunk_size = 100
    for i in range(0, len(issue_keys), chunk_size):
        chunk = issue_keys[i : i + chunk_size]
        print(f"{'[dry-run] ' if dry_run else ''}WONT_FIX {len(chunk)} issue(s) "
              f"(batch {i // chunk_size + 1})")
        if dry_run:
            for key in chunk[:5]:
                print(f"  e.g. {key}")
            if len(chunk) > 5:
                print(f"  ... and {len(chunk) - 5} more")
            continue
        _request(
            "POST",
            "/issues/bulk_change",
            token,
            data={
                "issues": ",".join(chunk),
                "do_transition": "wontfix",
                "comment": comment,
            },
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List matching open issues without resolving them",
    )
    parser.add_argument(
        "--reopen",
        metavar="RULES",
        help=(
            "Comma-separated rules to reopen instead of resolving, for when a "
            "recorded reason turns out to be wrong. A resolution nobody can "
            "defend is worse than an open finding."
        ),
    )
    args = parser.parse_args()
    token = os.environ.get("SONAR_TOKEN", "").strip()
    if not token:
        print(
            "SONAR_TOKEN is not set. Create a SonarCloud user token with "
            "Administer Issues on bckirkup_GutModelBacteriocins, then:\n"
            "  export SONAR_TOKEN=...\n"
            "  python3 scripts/sonar_wont_fix_debt.py",
            file=sys.stderr,
        )
        return 2

    if args.reopen:
        return reopen_rules(
            token,
            tuple(rule.strip() for rule in args.reopen.split(",")),
            dry_run=args.dry_run,
        )

    issues = list_open_debt_issues(token)
    print(f"Open issues matching accepted-debt rules: {len(issues)}")
    if not issues:
        print("Nothing to resolve.")
        return 0

    by_rule: dict[str, list[str]] = {}
    for issue in issues:
        by_rule.setdefault(issue["rule"], []).append(issue["key"])
    for rule, keys in sorted(
        by_rule.items(), key=lambda kv: (-len(kv[1]), kv[0])
    ):
        print(f"  {rule}: {len(keys)} — {DEBT_RULES[rule]}")

    # One transition per rule, so each resolution carries its own reason
    # instead of a single generic comment across unrelated rule families.
    for rule, keys in sorted(by_rule.items()):
        bulk_wont_fix(
            token,
            keys,
            f"{COMMENT_PREFIX}{rule}: {DEBT_RULES[rule]}",
            dry_run=args.dry_run,
        )
    if not args.dry_run:
        remaining = list_open_debt_issues(token)
        print(f"Remaining open accepted-debt issues: {len(remaining)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
