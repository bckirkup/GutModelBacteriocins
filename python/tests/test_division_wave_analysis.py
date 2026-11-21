"""Tests for division-wave-robust competitive strain-ratio estimators (issue #481)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from gut_ibm_tools import analysis


def test_ols_slope_recovers_trend_under_synchronized_division_waves() -> None:
    # Simulate 6 hours sampled every 10 minutes with 45-min synchronized division waves
    # producing +-0.20 decade single-step oscillations on top of a true +0.065 dec/6h
    # (+0.010833 dec/h) trend -- reproducing the b12=1e-5 aliasing landmine where the
    # final step alone reads +0.165 decades (~2.5x overstatement).
    times_s = np.arange(0.0, 21600.0 + 1.0, 600.0)
    true_slope_per_h = 0.065 / 6.0
    period_s = 2700.0
    phase_wave = 0.10 * np.sin(2.0 * np.pi * times_s / period_s + 0.5 * np.pi)
    # Shift baseline so t=0 starts near 0 and t=21600 peaks near +0.165
    log_ratio = true_slope_per_h * (times_s / 3600.0) + phase_wave

    n_comp = np.full_like(times_s, 1000.0)
    n_focal = (n_comp + 0.5) * np.power(10.0, log_ratio) - 0.5

    ols = analysis.log_ratio_ols_slope(
        times_s,
        n_focal,
        n_comp,
        pseudocount=0.5,
        doubling_time_s=period_s,
    )
    assert ols["slope_decades_per_h"] is not None
    assert math.isclose(
        float(ols["slope_decades_per_h"]),
        true_slope_per_h,
        abs_tol=0.004,
    )
    assert float(ols["residual_peak_to_trough_decades"]) > 0.15
    assert ols["wave_aliasing_warning"] is True

    # Over an integer multiple of wave periods (4 * 2700s = 10800s), trapezoid wave
    # averaging cancels the sinusoidal oscillation while single-frame endpoint biases high.
    wave = analysis.wave_averaged_log_ratio(
        times_s,
        n_focal,
        n_comp,
        pseudocount=0.5,
        tail_window_s=4.0 * period_s,
    )
    expected_mid_ratio = true_slope_per_h * ((21600.0 - 2.0 * period_s) / 3600.0)
    assert wave["trapezoid_mean_log10_ratio"] is not None
    assert math.isclose(
        float(wave["trapezoid_mean_log10_ratio"]),
        expected_mid_ratio,
        abs_tol=0.01,
    )
    assert float(wave["endpoint_log10_ratio"]) > float(wave["trapezoid_mean_log10_ratio"]) + 0.08
    assert float(wave["endpoint_bias_decades"]) > 0.08
    assert wave["wave_aliasing_warning"] is True


def test_ols_slope_clean_exponential_competition_no_warning() -> None:
    times_s = np.linspace(0.0, 21600.0, 37)
    true_slope_per_h = 0.33
    n_comp = np.full_like(times_s, 500.0)
    n_focal = (n_comp + 0.5) * np.power(10.0, true_slope_per_h * (times_s / 3600.0)) - 0.5

    res = analysis.log_ratio_ols_slope(times_s, n_focal, n_comp, doubling_time_s=2700.0)
    assert res["slope_decades_per_h"] is not None
    assert math.isclose(float(res["slope_decades_per_h"]), true_slope_per_h, rel_tol=1e-6)
    assert res["wave_aliasing_warning"] is False


def test_ols_slope_short_series_returns_warning() -> None:
    res = analysis.log_ratio_ols_slope([0.0], [10], [10])
    assert res["slope_decades_per_h"] is None
    assert res["wave_aliasing_warning"] is True

    with pytest.raises(ValueError, match="pseudocount"):
        analysis.log_ratio_ols_slope([0.0, 60.0], [10, 10], [10, 10], pseudocount=0.0)
