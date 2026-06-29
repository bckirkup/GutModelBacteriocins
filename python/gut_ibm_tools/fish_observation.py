"""
Simulated DNA-FISH and HCR-FISH observation models for experimental validation.

Per VADI §75 and issue #25, standard HiPR-FISH targeting immunity mRNA fails
because transcripts exist in single-digit copies per cell.  Multicopy plasmid
DNA-FISH and HCR-FISH amplification provide viable alternatives whose detection
limits, hybridization efficiency, and spatial resolution can be modeled from
GutIBM HDF5 output.

References:
    - VADI.md §75 (validation sensitivity)
    - docs/MECHANISMS.md (exclusion-radius validation approach)
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal

import numpy as np
from scipy.ndimage import gaussian_filter

from .analysis import monochromatic_patch_score
from .hdf5_reader import GutIBMData

DEFAULT_FISH_RNG_SEED = 0

# Typical multicopy targets in Enterobacteriaceae (literature order-of-magnitude).
DEFAULT_PLASMID_COPY_NUMBER = 25.0
DEFAULT_RRNA_OPERON_COPIES = 7.0
DEFAULT_IMMUNITY_MRNA_COPIES = 3.0

# Optical microscopy limits (m).
OPTICAL_PSF_SIGMA_M = 100e-9   # ~200 nm FWHM diffraction limit
SUPERRES_PSF_SIGMA_M = 25e-9   # ~50 nm effective resolution


class FishTechnique(str, Enum):
    """FISH modalities with different amplification and detection properties."""

    HIPR_FISH = "hipr_fish"
    DNA_FISH = "dna_fish"
    HCR_FISH = "hcr_fish"


@dataclass(frozen=True)
class FishProbe:
    """
    Probe definition for a molecular target.

    *target_kind* selects how per-agent copy numbers are inferred from HDF5:
      - ``plasmid`` — BI loci × plasmid copy number (multicopy DNA)
      - ``rrna`` — phylogroup rRNA operon copies (constant per cell)
      - ``immunity_mrna`` — basal immunity transcript copies (low)
    """

    name: str
    target_kind: Literal["plasmid", "rrna", "immunity_mrna"]
    technique: FishTechnique = FishTechnique.DNA_FISH
    hybridization_efficiency: float = 0.8
    copy_number: float = DEFAULT_PLASMID_COPY_NUMBER
    copy_number_std: float = 0.0
    hcr_amplification_factor: float = 1.0
    probe_brightness: float = 1.0
    detection_snr_threshold: float = 3.0

    def __post_init__(self) -> None:
        if not 0.0 <= self.hybridization_efficiency <= 1.0:
            raise ValueError("hybridization_efficiency must be in [0, 1]")
        if self.copy_number < 0:
            raise ValueError("copy_number must be non-negative")
        if self.hcr_amplification_factor < 1.0:
            raise ValueError("hcr_amplification_factor must be >= 1")


@dataclass(frozen=True)
class MicroscopyConfig:
    """Synthetic wide-field or super-resolution microscopy parameters."""

    image_size_px: tuple[int, int] = (256, 256)
    pixel_size_m: float = 0.5e-6
    psf_sigma_m: float = OPTICAL_PSF_SIGMA_M
    autofluorescence: float = 0.15
    read_noise_std: float = 0.05
    projection: Literal["xy", "xz", "yz"] = "xy"
    margin_m: float = 5e-6

    @property
    def resolution_m(self) -> float:
        """Approximate optical resolution (FWHM ≈ 2.355 σ)."""
        return 2.355 * self.psf_sigma_m

    @classmethod
    def optical(cls, **kwargs: object) -> "MicroscopyConfig":
        return cls(psf_sigma_m=OPTICAL_PSF_SIGMA_M, **kwargs)

    @classmethod
    def super_resolution(cls, **kwargs: object) -> "MicroscopyConfig":
        return cls(psf_sigma_m=SUPERRES_PSF_SIGMA_M, pixel_size_m=0.1e-6, **kwargs)


@dataclass
class FishObservationResult:
    """Per-agent and image-level outputs from a simulated FISH acquisition."""

    probe: FishProbe
    copy_numbers: np.ndarray
    raw_signal: np.ndarray
    snr: np.ndarray
    detected: np.ndarray
    detection_fraction: float
    mean_snr: float
    image: np.ndarray
    extent_m: tuple[float, float, float, float]


def _infer_copy_numbers(
    n_agents: int,
    probe: FishProbe,
    lineage: dict[str, np.ndarray],
    genome: dict[str, np.ndarray],
    rng: np.random.Generator,
) -> np.ndarray:
    """Map simulation state to per-agent target copy numbers."""
    if probe.target_kind == "plasmid":
        n_bi = lineage.get("num_bi_loci", np.zeros(n_agents, dtype=np.int32))
        base = np.maximum(n_bi.astype(np.float64), 0.0) * probe.copy_number
        # Agents without BI loci may still carry conjugative backbone plasmids.
        if "has_conjugative_plasmid" in genome:
            backbone = genome["has_conjugative_plasmid"].astype(np.float64)
            base = np.maximum(base, backbone * probe.copy_number * 0.5)
        copies = base
    elif probe.target_kind == "rrna":
        copies = np.full(n_agents, probe.copy_number)
    elif probe.target_kind == "immunity_mrna":
        copies = np.full(n_agents, probe.copy_number)
    else:
        raise ValueError(f"unknown target_kind: {probe.target_kind}")

    if probe.copy_number_std > 0:
        copies = copies + rng.normal(0.0, probe.copy_number_std, size=n_agents)
    return np.maximum(copies, 0.0)


def expected_probe_signal(
    copy_numbers: np.ndarray,
    probe: FishProbe,
) -> np.ndarray:
    """
    Expected fluorescence signal per agent before noise.

    Signal ∝ copy_number × hybridization_efficiency × amplification × brightness.
    """
    amp = probe.hcr_amplification_factor if probe.technique == FishTechnique.HCR_FISH else 1.0
    return (
        copy_numbers
        * probe.hybridization_efficiency
        * amp
        * probe.probe_brightness
    )


def compute_snr(
    signal: np.ndarray,
    *,
    background: float,
    read_noise_std: float,
) -> np.ndarray:
    """Per-agent signal-to-noise ratio against tissue autofluorescence."""
    noise = np.sqrt(np.maximum(background, 1e-12) + read_noise_std**2)
    return signal / noise


def detection_mask(snr: np.ndarray, threshold: float) -> np.ndarray:
    """Boolean mask of agents passing the SNR detection threshold."""
    return snr >= threshold


# Canonical probes for VADI validation scenarios.
PLASMID_DNA_FISH_PROBE = FishProbe(
    name="colicin_plasmid",
    target_kind="plasmid",
    technique=FishTechnique.DNA_FISH,
    hybridization_efficiency=0.85,
    copy_number=DEFAULT_PLASMID_COPY_NUMBER,
    copy_number_std=5.0,
    probe_brightness=1.0,
    detection_snr_threshold=3.0,
)

RRNA_PHYLOGROUP_PROBE = FishProbe(
    name="rrna_phylogroup",
    target_kind="rrna",
    technique=FishTechnique.DNA_FISH,
    hybridization_efficiency=0.9,
    copy_number=DEFAULT_RRNA_OPERON_COPIES,
    probe_brightness=1.2,
    detection_snr_threshold=3.0,
)

IMMUNITY_MRNA_HIPR_PROBE = FishProbe(
    name="immunity_mrna",
    target_kind="immunity_mrna",
    technique=FishTechnique.HIPR_FISH,
    hybridization_efficiency=0.6,
    copy_number=DEFAULT_IMMUNITY_MRNA_COPIES,
    copy_number_std=1.5,
    probe_brightness=0.5,
    detection_snr_threshold=3.0,
)

IMMUNITY_MRNA_HCR_PROBE = FishProbe(
    name="immunity_mrna_hcr",
    target_kind="immunity_mrna",
    technique=FishTechnique.HCR_FISH,
    hybridization_efficiency=0.7,
    copy_number=DEFAULT_IMMUNITY_MRNA_COPIES,
    copy_number_std=1.5,
    hcr_amplification_factor=500.0,
    probe_brightness=0.5,
    detection_snr_threshold=3.0,
)

DEFAULT_VALIDATION_PROBES: tuple[FishProbe, ...] = (
    PLASMID_DNA_FISH_PROBE,
    RRNA_PHYLOGROUP_PROBE,
    IMMUNITY_MRNA_HIPR_PROBE,
    IMMUNITY_MRNA_HCR_PROBE,
)


def _projection_axes(projection: str) -> tuple[int, int]:
    mapping = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}
    if projection not in mapping:
        raise ValueError(f"unknown projection: {projection}")
    return mapping[projection]


def render_synthetic_image(
    positions: np.ndarray,
    per_agent_signal: np.ndarray,
    detected: np.ndarray,
    config: MicroscopyConfig,
    *,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """
    Rasterize agent fluorescence into a 2D synthetic microscopy frame.

    Only *detected* agents contribute signal.  Each spot is blurred with a
    Gaussian PSF matching *config.psf_sigma_m*, then Poisson/Gaussian noise
    and autofluorescence are added.
    """
    if rng is None:
        rng = np.random.default_rng(DEFAULT_FISH_RNG_SEED)

    ax0, ax1 = _projection_axes(config.projection)
    coords = positions[:, [ax0, ax1]]
    if len(coords) == 0:
        h, w = config.image_size_px
        extent = (0.0, w * config.pixel_size_m, 0.0, h * config.pixel_size_m)
        return np.full((h, w), config.autofluorescence), extent

    margin = config.margin_m
    x_min, x_max = float(coords[:, 0].min() - margin), float(coords[:, 0].max() + margin)
    y_min, y_max = float(coords[:, 1].min() - margin), float(coords[:, 1].max() + margin)
    extent = (x_min, x_max, y_min, y_max)

    h, w = config.image_size_px
    image = np.full((h, w), config.autofluorescence, dtype=np.float64)

    x_span = max(x_max - x_min, 1e-12)
    y_span = max(y_max - y_min, 1e-12)

    for i in range(len(positions)):
        if not detected[i]:
            continue
        px = int((coords[i, 0] - x_min) / x_span * (w - 1))
        py = int((coords[i, 1] - y_min) / y_span * (h - 1))
        px = int(np.clip(px, 0, w - 1))
        py = int(np.clip(py, 0, h - 1))
        image[py, px] += per_agent_signal[i]

    sigma_px = config.psf_sigma_m / config.pixel_size_m
    if sigma_px > 0:
        image = gaussian_filter(image, sigma=sigma_px)

    # Shot noise on signal plus Gaussian read noise.
    signal_part = np.maximum(image - config.autofluorescence, 0.0)
    noisy_signal = rng.poisson(signal_part + 1e-9).astype(np.float64)
    image = config.autofluorescence + noisy_signal + rng.normal(
        0.0, config.read_noise_std, size=image.shape
    )
    return np.maximum(image, 0.0), extent


def simulate_fish_observation(
    positions: np.ndarray,
    probe: FishProbe,
    lineage: dict[str, np.ndarray],
    genome: dict[str, np.ndarray],
    microscopy: MicroscopyConfig | None = None,
    *,
    rng: np.random.Generator | None = None,
) -> FishObservationResult:
    """
    Full forward model: copy numbers → hybridization → SNR → detection → image.
    """
    if rng is None:
        rng = np.random.default_rng(DEFAULT_FISH_RNG_SEED)
    if microscopy is None:
        microscopy = MicroscopyConfig.optical(projection="xy")

    n_agents = len(positions)
    copies = _infer_copy_numbers(n_agents, probe, lineage, genome, rng)
    signal = expected_probe_signal(copies, probe)
    snr = compute_snr(
        signal,
        background=microscopy.autofluorescence,
        read_noise_std=microscopy.read_noise_std,
    )
    detected = detection_mask(snr, probe.detection_snr_threshold)

    image, extent = render_synthetic_image(
        positions, signal, detected, microscopy, rng=rng,
    )

    return FishObservationResult(
        probe=probe,
        copy_numbers=copies,
        raw_signal=signal,
        snr=snr,
        detected=detected,
        detection_fraction=float(np.mean(detected)) if n_agents else 0.0,
        mean_snr=float(np.mean(snr)) if n_agents else 0.0,
        image=image,
        extent_m=extent,
    )


def simulate_from_hdf5(
    data: GutIBMData,
    step: str,
    probe: FishProbe,
    microscopy: MicroscopyConfig | None = None,
    *,
    rng: np.random.Generator | None = None,
) -> FishObservationResult:
    """Run *simulate_fish_observation* on a GutIBM HDF5 step."""
    agents = data.get_agents(step)
    positions = np.column_stack([agents["x"], agents["y"], agents["z"]])
    lineage = data.get_lineage(step)
    genome = data.get_genome(step)
    return simulate_fish_observation(
        positions, probe, lineage, genome, microscopy, rng=rng,
    )


def compare_technique_detectability(
    positions: np.ndarray,
    lineage: dict[str, np.ndarray],
    genome: dict[str, np.ndarray],
    probes: tuple[FishProbe, ...] | None = None,
    *,
    rng: np.random.Generator | None = None,
) -> dict[str, dict[str, float]]:
    """
    Compare detection fractions across FISH modalities (issue #25 table).

    Returns per-probe summary metrics demonstrating that immunity mRNA is
    detectable with HCR-FISH and plasmid DNA-FISH but not standard HiPR-FISH.
    """
    if probes is None:
        probes = DEFAULT_VALIDATION_PROBES
    if rng is None:
        rng = np.random.default_rng(25)

    summary: dict[str, dict[str, float]] = {}
    for probe in probes:
        result = simulate_fish_observation(
            positions, probe, lineage, genome, rng=rng,
        )
        summary[probe.name] = {
            "technique": probe.technique.value,
            "target_kind": probe.target_kind,
            "detection_fraction": result.detection_fraction,
            "mean_snr": result.mean_snr,
            "mean_copy_number": float(np.mean(result.copy_numbers)),
            "hybridization_efficiency": probe.hybridization_efficiency,
            "amplification_factor": (
                probe.hcr_amplification_factor
                if probe.technique == FishTechnique.HCR_FISH
                else 1.0
            ),
        }
    return summary


def detected_spatial_clustering(
    positions: np.ndarray,
    types: np.ndarray,
    detected: np.ndarray,
    *,
    radius: float = 10e-6,
) -> dict[str, float]:
    """
    Monochromatic patch score computed only on FISH-detected agents.

    Bridges simulated microscopy back to exclusion-radius validation metrics.
    """
    mask = detected.astype(bool)
    n_detected = int(np.sum(mask))
    if n_detected < 2:
        return {
            "n_detected": float(n_detected),
            "detection_fraction": float(np.mean(detected)),
            "monochromatic_score_detected": 0.0,
        }

    score = monochromatic_patch_score(
        positions[mask], types[mask], radius=radius,
    )
    return {
        "n_detected": float(n_detected),
        "detection_fraction": float(np.mean(detected)),
        "monochromatic_score_detected": score,
    }


def validate_fish_observability(
    data: GutIBMData,
    step: str,
    *,
    probes: tuple[FishProbe, ...] | None = None,
    microscopy: MicroscopyConfig | None = None,
    rng: np.random.Generator | None = None,
) -> dict[str, float | dict[str, dict[str, float]]]:
    """
    Validation entry point for issue #25.

    Returns technique comparison plus spatial clustering on plasmid DNA-FISH
    detections (link to VADI exclusion-radius metrics).
    """
    agents = data.get_agents(step)
    positions = np.column_stack([agents["x"], agents["y"], agents["z"]])
    types = agents["type"]
    lineage = data.get_lineage(step)
    genome = data.get_genome(step)

    if rng is None:
        rng = np.random.default_rng(25)

    technique_summary = compare_technique_detectability(
        positions, lineage, genome, probes, rng=rng,
    )

    plasmid_result = simulate_from_hdf5(
        data, step, PLASMID_DNA_FISH_PROBE, microscopy, rng=rng,
    )
    spatial = detected_spatial_clustering(
        positions, types, plasmid_result.detected,
    )

    return {
        "technique_comparison": technique_summary,
        "plasmid_dna_fish_detection_fraction": plasmid_result.detection_fraction,
        "plasmid_dna_fish_mean_snr": plasmid_result.mean_snr,
        "immunity_hipr_detectable": float(
            technique_summary.get("immunity_mrna", {}).get("detection_fraction", 0.0) > 0.5
        ),
        "immunity_hcr_detectable": float(
            technique_summary.get("immunity_mrna_hcr", {}).get("detection_fraction", 0.0) > 0.5
        ),
        **spatial,
    }


def flatten_fish_metrics(
    result: dict[str, float | dict[str, dict[str, float]]],
) -> dict[str, float]:
    """Flatten *validate_fish_observability* output to scalar metrics for CI."""
    technique = result.get("technique_comparison", {})
    if not isinstance(technique, dict):
        technique = {}

    flat: dict[str, float] = {}
    for key, value in result.items():
        if key == "technique_comparison":
            continue
        flat[key] = float(value)  # type: ignore[arg-type]

    probe_keys = (
        ("immunity_mrna", "immunity_mrna_detection_fraction"),
        ("immunity_mrna_hcr", "immunity_mrna_hcr_detection_fraction"),
        ("colicin_plasmid", "colicin_plasmid_detection_fraction"),
        ("rrna_phylogroup", "rrna_phylogroup_detection_fraction"),
    )
    for probe_name, metric_name in probe_keys:
        probe_stats = technique.get(probe_name, {})
        if isinstance(probe_stats, dict) and "detection_fraction" in probe_stats:
            flat[metric_name] = float(probe_stats["detection_fraction"])

    return flat
