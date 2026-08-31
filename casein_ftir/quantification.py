"""
casein_ftir.quantification
==========================

Quantification :
    - Intégration des aires sous les bandes (avec correction de ligne
      de base locale)
    - Ratios diagnostiques A_amideII/A_amideI, A_phosphate/A_amideI, etc.
    - Application de la loi de Beer-Lambert avec étalonnage externe
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Optional, Tuple

from .io_module import Spectrum
from .database import CASEIN_BANDS, GALALITHE_MARKERS, IRBand
from .peak_detection import detect_bands, DetectedPeak


# ---------------------------------------------------------------------------
# Area integration
# ---------------------------------------------------------------------------

def integrate_band(spec: Spectrum, wn_min: float, wn_max: float,
                   baseline: str = "linear") -> float:
    """Integrated area under absorbance curve in [wn_min, wn_max].

    Parameters
    ----------
    baseline : {'linear', 'none'}
        'linear' : subtract a straight line between the two endpoints first.
        'none'   : raw integration.
    """
    mask = (spec.wavenumber >= wn_min) & (spec.wavenumber <= wn_max)
    if not np.any(mask):
        return 0.0
    x = spec.wavenumber[mask]
    y = spec.absorbance[mask]
    # Reorder ascending for integration
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if baseline == "linear" and len(x) > 1:
        y_base = np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
        y = y - y_base
        # clip negatives that may arise from overshoots
        y = np.clip(y, 0.0, None)

    return float(np.trapezoid(y, x))


def integrate_all_bands(spec: Spectrum, band_set: str = "casein"
                        ) -> Dict[str, float]:
    """Integrate every band of the chosen set."""
    if band_set == "casein":
        bands = CASEIN_BANDS
    elif band_set == "galalithe":
        bands = GALALITHE_MARKERS
    elif band_set == "all":
        bands = {**CASEIN_BANDS, **GALALITHE_MARKERS}
    else:
        raise ValueError(f"Unknown band_set: {band_set!r}")

    return {
        name: integrate_band(spec, band.window[0], band.window[1])
        for name, band in bands.items()
    }


# ---------------------------------------------------------------------------
# Diagnostic ratios
# ---------------------------------------------------------------------------

def diagnostic_ratios(spec: Spectrum) -> Dict[str, float]:
    """Compute key diagnostic ratios for casein.

    Returns
    -------
    dict with keys:
        - A_II / A_I              : amide II over amide I -- structural integrity
        - A_phosphate / A_I       : phosphorylation level
        - A_CH2(2926) / A_I       : aliphatic content vs protein backbone
        - A_amide_A / A_I         : N-H stretching strength
        - A_COO / A_I             : carboxylate content
    """
    areas = integrate_all_bands(spec, band_set="casein")
    a_I = areas.get("amide_I", 0.0)
    if a_I < 1e-12:
        return {k: np.nan for k in [
            "A_II_over_I", "A_phos_over_I", "A_CH2_over_I",
            "A_amideA_over_I", "A_COO_over_I",
        ]}

    return {
        "A_II_over_I":     areas.get("amide_II", 0.0) / a_I,
        "A_phos_over_I":   areas.get("phosphate", 0.0) / a_I,
        "A_CH2_over_I":    areas.get("CH2_asym", 0.0) / a_I,
        "A_amideA_over_I": areas.get("amide_A", 0.0) / a_I,
        "A_COO_over_I":    areas.get("COO_sym", 0.0) / a_I,
    }


def galalithe_indicators(spec: Spectrum) -> Dict[str, float]:
    """Indicators of casein -> galalithe polymerization.

    Returns
    -------
    dict with:
        - methylene_bridge_strength : area of ~2920 band
        - aminal_strength           : area of ~1370 band (N-CH2-N)
        - residual_formol           : area of ~1720 band (should be ~0)
        - crosslink_index           : composite indicator [0, +inf]
    """
    bridge = integrate_band(
        spec, *GALALITHE_MARKERS["methylene_bridge"].window
    )
    aminal = integrate_band(spec, *GALALITHE_MARKERS["N_CH2_N"].window)
    residual = integrate_band(
        spec, *GALALITHE_MARKERS["free_aldehyde_CO"].window
    )
    amide_I = integrate_band(spec, *CASEIN_BANDS["amide_I"].window)

    if amide_I < 1e-12:
        return {
            "methylene_bridge_strength": np.nan,
            "aminal_strength": np.nan,
            "residual_formol": np.nan,
            "crosslink_index": np.nan,
        }

    return {
        "methylene_bridge_strength": bridge / amide_I,
        "aminal_strength":           aminal / amide_I,
        "residual_formol":           residual / amide_I,
        "crosslink_index":           (bridge + aminal) / amide_I,
    }


# ---------------------------------------------------------------------------
# Beer-Lambert quantification
# ---------------------------------------------------------------------------

def beer_lambert_concentration(spec: Spectrum,
                                band_name: str,
                                molar_absorptivity: float,
                                path_length: float,
                                method: str = "peak") -> float:
    """Apply Beer-Lambert : A = epsilon * c * l  =>  c = A / (epsilon * l).

    Parameters
    ----------
    band_name : str
        Name of the band (must exist in CASEIN_BANDS).
    molar_absorptivity : float
        epsilon in L/(mol.cm). Must come from external calibration.
    path_length : float
        l in cm.
    method : {'peak', 'area'}
        'peak' : use the absorbance at the band maximum.
        'area' : use the integrated band area instead.

    Returns
    -------
    concentration : float
        In mol/L (or whichever unit matches epsilon).
    """
    if band_name not in CASEIN_BANDS:
        raise KeyError(f"Unknown casein band: {band_name!r}")
    band = CASEIN_BANDS[band_name]

    if method == "peak":
        peaks = detect_bands(spec, band_set="casein")
        p = peaks[band_name]
        if not p.found:
            return np.nan
        A = p.absorbance
    elif method == "area":
        A = integrate_band(spec, *band.window)
    else:
        raise ValueError(f"Unknown method: {method!r}")

    if molar_absorptivity * path_length < 1e-12:
        return np.nan
    return A / (molar_absorptivity * path_length)


def build_calibration_curve(spectra: list, concentrations: list,
                             band_name: str = "amide_II",
                             method: str = "area") -> Dict[str, float]:
    """Build a calibration curve from spectra of known concentrations.

    Linear regression : A = m * c + b

    Parameters
    ----------
    spectra : list of Spectrum
    concentrations : list of float (same length as spectra)
    band_name : str
        Band to use as quantification signal.
    method : {'peak', 'area'}

    Returns
    -------
    dict with keys 'slope', 'intercept', 'r_squared', 'band', 'method'
    """
    if len(spectra) != len(concentrations):
        raise ValueError("spectra and concentrations must have same length")
    if band_name not in CASEIN_BANDS:
        raise KeyError(f"Unknown band: {band_name!r}")
    band = CASEIN_BANDS[band_name]

    signals = []
    for s in spectra:
        if method == "area":
            signals.append(integrate_band(s, *band.window))
        else:
            peaks = detect_bands(s, band_set="casein")
            p = peaks[band_name]
            signals.append(p.absorbance if p.found else np.nan)

    c = np.array(concentrations, dtype=float)
    A = np.array(signals, dtype=float)
    valid = ~np.isnan(A)
    c, A = c[valid], A[valid]

    # Linear fit A = m c + b
    m, b = np.polyfit(c, A, 1)
    # R^2
    A_pred = m * c + b
    ss_res = np.sum((A - A_pred) ** 2)
    ss_tot = np.sum((A - A.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "slope": float(m),
        "intercept": float(b),
        "r_squared": float(r2),
        "band": band_name,
        "method": method,
        "n_points": int(valid.sum()),
    }


def apply_calibration(spec: Spectrum, calibration: Dict[str, float]) -> float:
    """Use a calibration curve to convert a spectrum's signal to concentration."""
    band_name = calibration["band"]
    method = calibration["method"]
    band = CASEIN_BANDS[band_name]
    if method == "area":
        A = integrate_band(spec, *band.window)
    else:
        peaks = detect_bands(spec, band_set="casein")
        p = peaks[band_name]
        A = p.absorbance if p.found else np.nan
    if np.isnan(A) or calibration["slope"] == 0:
        return float("nan")
    return (A - calibration["intercept"]) / calibration["slope"]
