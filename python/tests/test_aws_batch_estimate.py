"""Unit tests for AWS Batch cost estimates (no live AWS)."""

from __future__ import annotations

import pytest
from gut_ibm_tools.aws_batch_estimate import (
    PRICE_TABLE_USD_PER_HOUR,
    estimate_cost,
    format_estimate,
    lookup_price,
)


def test_lookup_price_spot_and_ondemand() -> None:
    spot = lookup_price("g4dn.xlarge", spot=True)
    od = lookup_price("g4dn.xlarge", spot=False)
    assert spot == pytest.approx(PRICE_TABLE_USD_PER_HOUR["g4dn.xlarge"]["spot"])
    assert od == pytest.approx(PRICE_TABLE_USD_PER_HOUR["g4dn.xlarge"]["on_demand"])
    assert spot < od


def test_estimate_campaign_scales_with_array() -> None:
    one = estimate_cost(instance_type="g5.2xlarge", wall_hours=10.0, array_size=1)
    twelve = estimate_cost(instance_type="g5.2xlarge", wall_hours=10.0, array_size=12)
    assert twelve.usd_campaign == pytest.approx(one.usd_per_run * 12)
    assert "usd_campaign" in format_estimate(twelve)


def test_interrupt_sensitivity_increases_cost() -> None:
    base = estimate_cost(instance_type="g5.2xlarge", wall_hours=24.0, spot=True)
    with_int = estimate_cost(
        instance_type="g5.2xlarge",
        wall_hours=24.0,
        spot=True,
        interrupt_every_hours=6.0,
        checkpoint_overhead_hours=0.25,
    )
    assert with_int.usd_per_run > base.usd_per_run


def test_unknown_instance_raises() -> None:
    with pytest.raises(KeyError):
        lookup_price("t2.micro", spot=True)
