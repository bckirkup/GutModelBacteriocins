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
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

PROJECT = "bckirkup_GutModelBacteriocins"
BASE = "https://sonarcloud.io/api"

# Accepted complexity / architecture debt (Batch B). Do not include rules
# cleared by mechanical code fixes (Batch A).
DEBT_RULES = (
    "cpp:S134",
    "cpp:S107",
    "cpp:S6004",
    "cpp:S3776",
    "python:S3776",
    "cpp:S995",
    "cpp:S5008",
    "cpp:S1820",
    "cpp:S1448",
    "cpp:S3656",
    "cpp:S924",
    "cpp:S7034",
)

COMMENT = (
    "Accepted GutIBM technical debt — see docs/SONARQUBE_PLAN.md. "
    "Complexity/architecture smells in QSSA/FMM/Simulation/Fix; "
    "multicriteria already suppresses for scanner runs."
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
            import json

            return json.loads(raw)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")
        raise SystemExit(f"HTTP {exc.code} on {path}: {detail}") from exc


def list_open_debt_issues(token: str) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    page = 1
    while True:
        payload = _request(
            "GET",
            "/issues/search",
            token,
            params={
                "componentKeys": PROJECT,
                "resolved": "false",
                "rules": ",".join(DEBT_RULES),
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


def bulk_wont_fix(token: str, issue_keys: list[str], dry_run: bool) -> None:
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
                "comment": COMMENT,
            },
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List matching open issues without resolving them",
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

    issues = list_open_debt_issues(token)
    print(f"Open debt issues matching Batch B rules: {len(issues)}")
    if not issues:
        print("Nothing to resolve.")
        return 0

    by_rule: dict[str, int] = {}
    for issue in issues:
        by_rule[issue["rule"]] = by_rule.get(issue["rule"], 0) + 1
    for rule, count in sorted(by_rule.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {rule}: {count}")

    bulk_wont_fix(token, [i["key"] for i in issues], dry_run=args.dry_run)
    if not args.dry_run:
        remaining = list_open_debt_issues(token)
        print(f"Remaining open debt issues: {len(remaining)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
