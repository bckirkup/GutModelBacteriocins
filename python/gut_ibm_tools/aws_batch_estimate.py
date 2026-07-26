"""Rough AWS Batch $/run estimates for GutIBM GPU jobs.

Uses a checked-in us-east-1 price table. Optional live Spot price lookup when
credentials exist (injected / ``aws ec2 describe-spot-price-history``).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass

AWS = "aws"
DEFAULT_REGION = "us-east-1"

# Approximate us-east-1 list prices (USD/hour). Spot midpoints are illustrative;
# refresh from AWS Pricing / Spot history when planning a campaign.
PRICE_TABLE_USD_PER_HOUR: dict[str, dict[str, float]] = {
    "g4dn.xlarge": {"on_demand": 0.526, "spot": 0.18},
    "g4dn.2xlarge": {"on_demand": 0.752, "spot": 0.26},
    "g5.2xlarge": {"on_demand": 1.212, "spot": 0.42},
    "g6.2xlarge": {"on_demand": 0.978, "spot": 0.35},
}

AwsRunner = Callable[[list[str]], str]


@dataclass(frozen=True)
class CostEstimate:
    instance_type: str
    region: str
    spot: bool
    usd_per_hour: float
    wall_hours: float
    array_size: int
    interrupt_every_hours: float | None
    checkpoint_overhead_hours: float
    usd_per_run: float
    usd_campaign: float
    notes: str


def lookup_price(
    instance_type: str,
    *,
    spot: bool,
    region: str = DEFAULT_REGION,
    price_table: dict[str, dict[str, float]] | None = None,
) -> float:
    table = price_table or PRICE_TABLE_USD_PER_HOUR
    if instance_type not in table:
        known = ", ".join(sorted(table))
        raise KeyError(f"unknown instance type {instance_type!r}; known: {known}")
    key = "spot" if spot else "on_demand"
    return float(table[instance_type][key])


def fetch_live_spot_price(
    instance_type: str,
    *,
    region: str = DEFAULT_REGION,
    aws_run: AwsRunner | None = None,
) -> float | None:
    """Return the latest Spot price or None if the call fails / is empty."""
    runner = aws_run or (
        lambda argv: subprocess.run(
            argv, capture_output=True, text=True, check=True
        ).stdout
    )
    try:
        raw = runner(
            [
                AWS,
                "ec2",
                "describe-spot-price-history",
                "--region",
                region,
                "--instance-types",
                instance_type,
                "--product-descriptions",
                "Linux/UNIX",
                "--max-items",
                "1",
                "--output",
                "json",
            ]
        )
        hist = json.loads(raw).get("SpotPriceHistory") or []
        if not hist:
            return None
        return float(hist[0]["SpotPrice"])
    except (subprocess.CalledProcessError, KeyError, TypeError, ValueError, OSError):
        return None


def estimate_cost(
    *,
    instance_type: str,
    wall_hours: float,
    array_size: int = 1,
    spot: bool = True,
    region: str = DEFAULT_REGION,
    interrupt_every_hours: float | None = None,
    checkpoint_overhead_hours: float = 0.1,
    usd_per_hour: float | None = None,
    price_table: dict[str, dict[str, float]] | None = None,
) -> CostEstimate:
    """Estimate $/run and $/campaign with optional Spot interrupt sensitivity."""
    if wall_hours < 0 or array_size < 1:
        raise ValueError("wall_hours must be >= 0 and array_size >= 1")
    rate = (
        usd_per_hour
        if usd_per_hour is not None
        else lookup_price(instance_type, spot=spot, region=region, price_table=price_table)
    )
    effective_hours = wall_hours
    notes = "table price"
    if interrupt_every_hours is not None and interrupt_every_hours > 0:
        # Rough: each reclaim wastes up to one checkpoint interval of work plus
        # restart overhead; count expected interrupts over the wall clock.
        n_interrupt = wall_hours / interrupt_every_hours
        effective_hours = wall_hours + n_interrupt * checkpoint_overhead_hours
        notes = (
            f"table price + {n_interrupt:.1f} interrupts "
            f"× {checkpoint_overhead_hours}h overhead"
        )
    usd_per_run = rate * effective_hours
    return CostEstimate(
        instance_type=instance_type,
        region=region,
        spot=spot,
        usd_per_hour=rate,
        wall_hours=wall_hours,
        array_size=array_size,
        interrupt_every_hours=interrupt_every_hours,
        checkpoint_overhead_hours=checkpoint_overhead_hours,
        usd_per_run=usd_per_run,
        usd_campaign=usd_per_run * array_size,
        notes=notes,
    )


def format_estimate(est: CostEstimate) -> str:
    market = "Spot" if est.spot else "On-Demand"
    return "\n".join(
        [
            f"instance:        {est.instance_type} ({market}) @ {est.region}",
            f"usd_per_hour:    {est.usd_per_hour:.4f}",
            f"wall_hours:      {est.wall_hours}",
            f"array_size:      {est.array_size}",
            f"usd_per_run:     {est.usd_per_run:.2f}",
            f"usd_campaign:    {est.usd_campaign:.2f}",
            f"notes:           {est.notes}",
        ]
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Rough GutIBM AWS Batch GPU cost estimate (us-east-1 table).",
    )
    parser.add_argument(
        "--instance-type",
        default="g5.2xlarge",
        help="EC2 instance type (must be in the checked-in price table)",
    )
    parser.add_argument(
        "--wall-hours",
        type=float,
        required=True,
        help="Expected wall-clock hours per simulation run",
    )
    parser.add_argument("--array-size", type=int, default=1)
    parser.add_argument(
        "--on-demand",
        action="store_true",
        help="Use On-Demand table price instead of Spot",
    )
    parser.add_argument("--region", default=DEFAULT_REGION)
    parser.add_argument(
        "--interrupt-every-hours",
        type=float,
        default=None,
        help="Sensitivity: assume Spot reclaim every N wall hours",
    )
    parser.add_argument(
        "--checkpoint-overhead-hours",
        type=float,
        default=0.1,
        help="Extra hours charged per interrupt (restart / lost work)",
    )
    parser.add_argument(
        "--live-spot",
        action="store_true",
        help="Try aws ec2 describe-spot-price-history (falls back to table)",
    )
    parser.add_argument(
        "--list-prices",
        action="store_true",
        help="Print the checked-in price table and exit",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.list_prices:
        print(json.dumps(PRICE_TABLE_USD_PER_HOUR, indent=2))
        return 0
    live: float | None = None
    if args.live_spot and not args.on_demand:
        live = fetch_live_spot_price(args.instance_type, region=args.region)
    try:
        est = estimate_cost(
            instance_type=args.instance_type,
            wall_hours=args.wall_hours,
            array_size=args.array_size,
            spot=not args.on_demand,
            region=args.region,
            interrupt_every_hours=args.interrupt_every_hours,
            checkpoint_overhead_hours=args.checkpoint_overhead_hours,
            usd_per_hour=live,
        )
    except (KeyError, ValueError) as exc:
        print(f"aws batch estimate error: {exc}", file=sys.stderr)
        return 2
    if live is not None:
        print(f"(live Spot price used: {live:.4f} USD/h)")
    print(format_estimate(est))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
