"""Tests for HCR-FISH / DNA-FISH observation models (issue #25)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from gut_ibm_tools import GutIBMData, fish_observation
from gut_ibm_tools.fish_observation import (
    FishProbe,
    FishTechnique,
    MicroscopyConfig,
    PLASMID_DNA_FISH_PROBE,
    compute_snr,
    detection_mask,
    expected_probe_signal,
    simulate_fish_observation,
    validate_fish_observability,
)


def _synthetic_agents(
    n_per_type: int = 20,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray]:
    rng = np.random.default_rng(7)
    n = n_per_type * 2
    positions = np.vstack([
        rng.uniform(0, 30e-6, size=(n_per_type, 3)),
        rng.uniform(70e-6, 100e-6, size=(n_per_type, 3)),
    ])
    types = np.array([1] * n_per_type + [2] * n_per_type, dtype=np.int32)
    lineage = {
        "num_bi_loci": np.where(types == 1, 2, 1).astype(np.int32),
        "btuB_expression": rng.uniform(0.3, 1.0, n),
    }
    genome = {
        "has_conjugative_plasmid": np.ones(n, dtype=np.int32),
    }
    return positions, lineage, genome, types


def test_expected_probe_signal_scales_with_copies() -> None:
    copies = np.array([0.0, 10.0, 50.0])
    probe = FishProbe(
        name="test",
        target_kind="plasmid",
        hybridization_efficiency=0.5,
        copy_number=1.0,
        probe_brightness=2.0,
    )
    signal = expected_probe_signal(copies, probe)
    np.testing.assert_allclose(signal, copies * 0.5 * 2.0)


def test_hcr_amplification_boosts_low_copy_signal() -> None:
    copies = np.full(5, 3.0)
    hipr = FishProbe(
        name="hipr",
        target_kind="immunity_mrna",
        technique=FishTechnique.HIPR_FISH,
        copy_number=3.0,
        hybridization_efficiency=0.6,
    )
    hcr = FishProbe(
        name="hcr",
        target_kind="immunity_mrna",
        technique=FishTechnique.HCR_FISH,
        copy_number=3.0,
        hybridization_efficiency=0.6,
        hcr_amplification_factor=500.0,
    )
    assert np.mean(expected_probe_signal(copies, hcr)) > np.mean(
        expected_probe_signal(copies, hipr),
    ) * 100


def test_immunity_mrna_hipr_fails_hcr_succeeds() -> None:
    positions, lineage, genome, _ = _synthetic_agents()
    rng = np.random.default_rng(25)
    summary = fish_observation.compare_technique_detectability(
        positions, lineage, genome, rng=rng,
    )

    hipr_frac = summary["immunity_mrna"]["detection_fraction"]
    hcr_frac = summary["immunity_mrna_hcr"]["detection_fraction"]
    plasmid_frac = summary["colicin_plasmid"]["detection_fraction"]

    assert hipr_frac < 0.5
    assert hcr_frac > 0.9
    assert plasmid_frac > 0.9


def test_plasmid_dna_fish_detects_multicopy_targets() -> None:
    positions, lineage, genome, _ = _synthetic_agents()
    result = simulate_fish_observation(
        positions,
        PLASMID_DNA_FISH_PROBE,
        lineage,
        genome,
        MicroscopyConfig.optical(image_size_px=(128, 128)),
        rng=np.random.default_rng(0),
    )

    assert result.detection_fraction > 0.9
    assert result.image.shape == (128, 128)
    assert result.image.max() > result.image.min()
    assert result.mean_snr > PLASMID_DNA_FISH_PROBE.detection_snr_threshold


def test_render_synthetic_image_resolution_limits() -> None:
    positions = np.array([[0.0, 0.0, 0.0], [0.2e-6, 0.0, 0.0]])
    signal = np.array([10.0, 10.0])
    detected = np.array([True, True])

    optical = MicroscopyConfig.optical(
        image_size_px=(64, 64),
        pixel_size_m=0.05e-6,
        margin_m=1e-6,
    )
    superres = MicroscopyConfig.super_resolution(
        image_size_px=(64, 64),
        margin_m=1e-6,
    )

    img_opt, _ = fish_observation.render_synthetic_image(
        positions, signal, detected, optical, rng=np.random.default_rng(1),
    )
    img_sr, _ = fish_observation.render_synthetic_image(
        positions, signal, detected, superres, rng=np.random.default_rng(1),
    )

    assert optical.resolution_m > superres.resolution_m
    # Super-resolution PSF is narrower → peak should be sharper (higher max).
    assert img_sr.max() >= img_opt.max() * 0.5


def test_detected_spatial_clustering_segregated_patches() -> None:
    positions, _, _, types = _synthetic_agents()
    detected = np.ones(len(positions), dtype=bool)

    metrics = fish_observation.detected_spatial_clustering(
        positions, types, detected, radius=15e-6,
    )
    assert metrics["monochromatic_score_detected"] > 0.7


def test_validate_fish_observability_from_hdf5(sample_hdf5: Path) -> None:
    with GutIBMData(sample_hdf5) as data:
        result = validate_fish_observability(data, "step_000000")

    assert "technique_comparison" in result
    assert result["plasmid_dna_fish_detection_fraction"] > 0.5
    assert result["immunity_hipr_detectable"] == pytest.approx(0.0)
    assert result["immunity_hcr_detectable"] == pytest.approx(1.0)
    assert result["monochromatic_score_detected"] > 0.5


def test_compute_snr_and_detection_mask() -> None:
    signal = np.array([0.0, 3.0, 9.0])
    snr = compute_snr(signal, background=1.0, read_noise_std=0.1)
    detected = detection_mask(snr, threshold=3.0)
    assert detected.tolist() == [False, False, True]


def test_flatten_fish_metrics_roundtrip() -> None:
    from gut_ibm_tools.fish_observation import flatten_fish_metrics

    raw = {
        "technique_comparison": {
            "immunity_mrna": {"detection_fraction": 0.1},
            "colicin_plasmid": {"detection_fraction": 0.9},
        },
        "plasmid_dna_fish_detection_fraction": 0.9,
        "immunity_hipr_detectable": 0.0,
    }
    flat = flatten_fish_metrics(raw)
    assert flat["immunity_mrna_detection_fraction"] == pytest.approx(0.1)
    assert flat["colicin_plasmid_detection_fraction"] == pytest.approx(0.9)
    assert flat["plasmid_dna_fish_detection_fraction"] == pytest.approx(0.9)


def test_fish_probe_validation_errors() -> None:
    with pytest.raises(ValueError, match="hybridization_efficiency"):
        FishProbe(
            name="bad",
            target_kind="rrna",
            hybridization_efficiency=1.5,
        )
