"""
casein_ftir.simulate
====================

Simulation paramétrique de spectres FTIR de caséine et de galalithe.

Méthode : somme de gaussiennes positionnées aux fréquences caractéristiques
issues de la littérature (cf. casein_ftir.database). Les paramètres sont
ajustables : positions, intensités, largeurs, bruit.

PROPOSITION — Modèle pédagogique par sommes de bandes. Les déplacements,
largeurs, intensités et paramètres de réticulation sont des hypothèses.
Ce modèle n'est ni un calcul ab initio, ni un étalon expérimental.
À VALIDER — Toute interprétation des paramètres comme propriétés du matériau.
"""

from __future__ import annotations
import numpy as np
from typing import Optional, Dict, Iterable
from .io_module import Spectrum
from .database import CASEIN_BANDS, AMIDE_I_SUBBANDS, GALALITHE_MARKERS, IRBand


# ---------------------------------------------------------------------------
# Line shapes
# ---------------------------------------------------------------------------

def gaussian(x: np.ndarray, center: float, amplitude: float,
             fwhm: float) -> np.ndarray:
    """Gaussian line shape."""
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def lorentzian(x: np.ndarray, center: float, amplitude: float,
               fwhm: float) -> np.ndarray:
    """Lorentzian line shape."""
    gamma = fwhm / 2.0
    return amplitude * gamma ** 2 / ((x - center) ** 2 + gamma ** 2)


def voigt(x: np.ndarray, center: float, amplitude: float,
          fwhm_g: float, fwhm_l: float) -> np.ndarray:
    """Pseudo-Voigt (50/50 mix). For real Voigt use lmfit."""
    return 0.5 * gaussian(x, center, amplitude, fwhm_g) + \
           0.5 * lorentzian(x, center, amplitude, fwhm_l)


# ---------------------------------------------------------------------------
# Species-specific tweaks (small literature-based shifts)
# ---------------------------------------------------------------------------

SPECIES_SHIFTS: Dict[str, Dict[str, float]] = {
    # Cow casein -- reference, zero shifts
    "bovine": {},
    # Goat: alpha-s1 reduced, slight shifts due to different secondary
    # structure proportions. Values approximated from literature.
    "caprine": {"amide_I": -2.0, "phosphate": +5.0},
    # Sheep: composition closer to cow but higher protein content.
    "ovine":   {"amide_I": -1.0, "amide_II": +2.0},
    # Buffalo, human, etc. -- placeholder for future extension
}

DEFAULT_FWHM = {
    # Reasonable FWHM (cm^-1) for ATR-FTIR on solid casein
    "amide_A": 200.0,
    "CH3_asym": 20.0,
    "CH2_asym": 20.0,
    "CH_sym": 20.0,
    "amide_I": 60.0,
    "amide_II": 50.0,
    "CH2_bend": 25.0,
    "COO_sym": 35.0,
    "amide_III": 60.0,
    "phosphate": 80.0,
    "C_O_skeletal": 90.0,
    "methylene_bridge": 25.0,
    "free_aldehyde_CO": 30.0,
    "N_CH2_N": 40.0,
}


# ---------------------------------------------------------------------------
# Casein simulator
# ---------------------------------------------------------------------------

def simulate_casein(wn_min: float = 600.0,
                    wn_max: float = 4000.0,
                    resolution: float = 2.0,
                    species: str = "bovine",
                    secondary_structure: Optional[Dict[str, float]] = None,
                    noise_level: float = 0.0,
                    baseline_drift: float = 0.0) -> Spectrum:
    """Generate a synthetic FTIR spectrum of casein.

    Parameters
    ----------
    wn_min, wn_max : float
        Spectral range in cm^-1.
    resolution : float
        Point spacing in cm^-1 (typical FTIR: 2 or 4).
    species : {'bovine', 'caprine', 'ovine'}
        Applies illustrative, unvalidated peak shifts.
    secondary_structure : dict, optional
        Fractions of {'alpha_helix', 'beta_sheet', 'random_coil',
        'beta_turn', 'beta_aggregate', 'beta_antiparallel'} summing to 1.0.
        If given, the amide I band is split into the corresponding sub-bands
        with these relative weights. Default : casein-like (mostly random
        coil, see Barth 2007 and D'Incecco 2025).
    noise_level : float
        Gaussian white noise stddev added to spectrum (in absorbance units).
    baseline_drift : float
        Amplitude of a low-frequency sinusoidal baseline drift, to simulate
        instrumental artifacts.

    Returns
    -------
    Spectrum
    """
    if not all(np.isfinite(v) for v in (wn_min, wn_max, resolution, noise_level, baseline_drift)):
        raise ValueError("Simulation parameters must be finite")
    if wn_min >= wn_max or resolution <= 0 or noise_level < 0 or baseline_drift < 0:
        raise ValueError("Invalid range, resolution, noise or baseline parameter")
    if species not in SPECIES_SHIFTS:
        raise ValueError(f"Unknown species: {species}")
    # Build descending wavenumber axis (IR convention)
    wn = np.arange(wn_max, wn_min - resolution, -resolution)
    abs_ = np.zeros_like(wn)

    shifts = SPECIES_SHIFTS.get(species, {})

    # --- Default secondary structure for native casein ----------------------
    if secondary_structure is None:
        # Approx. from Barth 2007 + D'Incecco 2025 (random-coil-dominated):
        secondary_structure = {
            "random_coil":       0.40,
            "alpha_helix":       0.20,
            "beta_sheet":        0.15,
            "beta_turn":         0.10,
            "beta_aggregate":    0.10,
            "beta_antiparallel": 0.05,
        }
    # Total amide I area = 1.0 (reference)
    amide_I_total = 1.0
    for sub_name, frac in secondary_structure.items():
        if sub_name not in AMIDE_I_SUBBANDS:
            continue
        sub = AMIDE_I_SUBBANDS[sub_name]
        center = sub.center + shifts.get("amide_I", 0.0)
        amplitude = amide_I_total * frac
        # narrower FWHM for sub-bands than the overall amide I envelope
        abs_ += gaussian(wn, center, amplitude, fwhm=25.0)

    # --- Other casein bands -------------------------------------------------
    for name, band in CASEIN_BANDS.items():
        if name == "amide_I":
            continue  # already handled via sub-bands
        center = band.center + shifts.get(name, 0.0)
        amplitude = band.rel_intensity_ref
        fwhm = DEFAULT_FWHM.get(name, 30.0)
        abs_ += gaussian(wn, center, amplitude, fwhm)

    # --- Optional artefacts -------------------------------------------------
    if noise_level > 0:
        rng = np.random.default_rng(seed=42)
        abs_ = abs_ + rng.normal(0.0, noise_level, size=abs_.shape)

    if baseline_drift > 0:
        # gentle sinusoidal drift over full range
        period = (wn_max - wn_min)
        abs_ = abs_ + baseline_drift * np.sin(2 * np.pi * (wn - wn_min) / period)

    return Spectrum(
        wavenumber=wn,
        absorbance=abs_,
        name=f"simulated_casein_{species}",
        metadata={
            "type": "simulated",
            "species": species,
            "secondary_structure": secondary_structure,
            "noise_level": noise_level,
            "baseline_drift": baseline_drift,
        },
    )


# ---------------------------------------------------------------------------
# Galalithe simulator (casein + crosslinks)
# ---------------------------------------------------------------------------

def simulate_galalithe(wn_min: float = 600.0,
                       wn_max: float = 4000.0,
                       resolution: float = 2.0,
                       crosslink_degree: float = 0.5,
                       residual_formol: float = 0.05,
                       species: str = "bovine",
                       noise_level: float = 0.0) -> Spectrum:
    """Simulate a polymerized casein-formaldehyde spectrum (galalithe).

    Built as a base casein spectrum + crosslinking markers proportional
    to ``crosslink_degree`` in [0, 1].

    Parameters
    ----------
    crosslink_degree : float in [0, 1]
        Synthetic marker amplitude only; not a measured conversion.
    residual_formol : float in [0, 1]
        Synthetic carbonyl amplitude only; not a formaldehyde concentration.
    """
    if not all(np.isfinite(v) and 0 <= v <= 1 for v in (crosslink_degree, residual_formol)):
        raise ValueError("Synthetic marker parameters must lie in [0, 1]")
    # Start from a casein spectrum -- but shift its secondary structure
    # towards more beta-aggregate (typical after crosslinking)
    crosslinked_struct = {
        "random_coil":       0.40 - 0.15 * crosslink_degree,
        "alpha_helix":       0.20 - 0.05 * crosslink_degree,
        "beta_sheet":        0.15 + 0.05 * crosslink_degree,
        "beta_turn":         0.10,
        "beta_aggregate":    0.10 + 0.10 * crosslink_degree,
        "beta_antiparallel": 0.05 + 0.05 * crosslink_degree,
    }
    base = simulate_casein(
        wn_min=wn_min, wn_max=wn_max, resolution=resolution,
        species=species, secondary_structure=crosslinked_struct,
        noise_level=0.0, baseline_drift=0.0,
    )

    # Add crosslink markers
    wn = base.wavenumber
    abs_ = base.absorbance.copy()

    # 1) Intensification du pic CH2 a ~2920 cm-1 (ponts methylene)
    bridge = GALALITHE_MARKERS["methylene_bridge"]
    abs_ += gaussian(wn, bridge.center,
                     amplitude=bridge.rel_intensity_ref * crosslink_degree,
                     fwhm=DEFAULT_FWHM["methylene_bridge"])

    # 2) Aminals N-CH2-N a ~1370 cm-1
    aminal = GALALITHE_MARKERS["N_CH2_N"]
    abs_ += gaussian(wn, aminal.center,
                     amplitude=aminal.rel_intensity_ref * crosslink_degree,
                     fwhm=DEFAULT_FWHM["N_CH2_N"])

    # 3) Residual free formaldehyde C=O (should ideally be near zero)
    if residual_formol > 0:
        free = GALALITHE_MARKERS["free_aldehyde_CO"]
        abs_ += gaussian(wn, free.center,
                         amplitude=residual_formol * 0.5,
                         fwhm=DEFAULT_FWHM["free_aldehyde_CO"])

    if noise_level > 0:
        rng = np.random.default_rng(seed=43)
        abs_ = abs_ + rng.normal(0.0, noise_level, size=abs_.shape)

    return Spectrum(
        wavenumber=wn,
        absorbance=abs_,
        name=f"simulated_galalithe_x{crosslink_degree:.2f}",
        metadata={
            "type": "simulated",
            "subtype": "galalithe",
            "species": species,
            "crosslink_degree": crosslink_degree,
            "residual_formol": residual_formol,
        },
    )
