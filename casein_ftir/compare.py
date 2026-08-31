"""
casein_ftir.compare
===================

Comparaison de deux spectres FTIR (typiquement avant/après polymérisation
caséine -> galalithe).

Fonctions :
    - resample_to_common_axis : ré-échantillonnage sur axe commun
    - difference_spectrum     : spectre de différence (apres - avant)
    - highlight_changes       : identification des bandes apparues/disparues
    - similarity_metric       : score de similarité global (corrélation)
"""

from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple

from .io_module import Spectrum
from .database import CASEIN_BANDS, GALALITHE_MARKERS
from .quantification import integrate_band


def resample_to_common_axis(s1: Spectrum, s2: Spectrum,
                             resolution: float = 2.0
                             ) -> Tuple[Spectrum, Spectrum]:
    """Re-sample both spectra on a common descending wavenumber axis."""
    wn_min = max(s1.wavenumber.min(), s2.wavenumber.min())
    wn_max = min(s1.wavenumber.max(), s2.wavenumber.max())
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("resolution must be finite and positive")
    if wn_max - wn_min < resolution:
        raise ValueError("Spectra need an overlapping range with at least two resampled points")
    n_points = int(np.floor((wn_max - wn_min) / resolution)) + 1
    wn_new = wn_max - resolution * np.arange(n_points)

    def _interp(s, x_new):
        # ascending for np.interp
        x_asc = s.wavenumber[::-1]
        y_asc = s.absorbance[::-1]
        return np.interp(x_new[::-1], x_asc, y_asc)[::-1]

    a1 = _interp(s1, wn_new)
    a2 = _interp(s2, wn_new)
    out1 = Spectrum(wn_new, a1, name=s1.name)
    out2 = Spectrum(wn_new, a2, name=s2.name)
    return out1, out2


def difference_spectrum(before: Spectrum, after: Spectrum,
                        normalize: bool = True) -> Spectrum:
    """Compute the absorbance difference ``after - before``.

    Positive peaks : bands appeared or grew (e.g., methylene bridges
    in galalithe). Negative peaks : bands shrank or disappeared (e.g.,
    free amine consumed by formaldehyde).
    """
    b, a = resample_to_common_axis(before, after)
    if normalize:
        # Normalise both on their amide I peak before subtraction.
        from .preprocess import normalize_to_band
        b = normalize_to_band(b, center=1648.0, halfwidth=30.0)
        a = normalize_to_band(a, center=1648.0, halfwidth=30.0)

    diff = Spectrum(
        wavenumber=a.wavenumber.copy(),
        absorbance=a.absorbance - b.absorbance,
        name=f"({a.name})-({b.name})",
        metadata={"type": "difference"},
    )
    return diff


def highlight_changes(before: Spectrum, after: Spectrum,
                      threshold: float = 0.05) -> Dict[str, Dict]:
    """For each known casein/galalithe band, report area change.

    Parameters
    ----------
    threshold : float
        Minimum relative change (|Delta/A_before|) to flag a band as
        "changed". 0.05 = 5%.

    Returns
    -------
    dict { band_name : {'A_before', 'A_after', 'delta', 'rel_change', 'flag'} }
    """
    all_bands = {**CASEIN_BANDS, **GALALITHE_MARKERS}
    out = {}
    for name, band in all_bands.items():
        A_b = integrate_band(before, *band.window)
        A_a = integrate_band(after,  *band.window)
        delta = A_a - A_b
        if A_b > 1e-9:
            rel = delta / A_b
        else:
            rel = float("inf") if A_a > 1e-9 else 0.0

        if rel == float("inf"):
            flag = "APPEARED"
        elif abs(rel) < threshold:
            flag = "stable"
        elif rel > 0:
            flag = "increased"
        else:
            flag = "decreased"

        out[name] = {
            "A_before":   A_b,
            "A_after":    A_a,
            "delta":      delta,
            "rel_change": rel,
            "flag":       flag,
        }
    return out


def similarity_metric(s1: Spectrum, s2: Spectrum) -> Dict[str, float]:
    """Several similarity scores between two spectra.

    Returns
    -------
    dict with:
        - pearson_corr : Pearson correlation coefficient
        - cosine       : cosine similarity (vector angle)
        - rmse         : root mean square error
        - spectral_overlap : 1 - normalised area of |A-B|
    """
    a, b = resample_to_common_axis(s1, s2)
    x = a.absorbance
    y = b.absorbance

    # Pearson
    if x.std() > 0 and y.std() > 0:
        pearson = float(np.corrcoef(x, y)[0, 1])
    else:
        pearson = float("nan")

    # Cosine
    nx = float(np.linalg.norm(x))
    ny = float(np.linalg.norm(y))
    cosine = float(np.dot(x, y) / (nx * ny)) if nx * ny > 0 else float("nan")

    # RMSE
    rmse = float(np.sqrt(np.mean((x - y) ** 2)))

    # Spectral overlap
    diff_area = float(np.trapezoid(np.abs(x - y), -a.wavenumber))
    tot_area  = float(np.trapezoid(0.5 * (np.abs(x) + np.abs(y)), -a.wavenumber))
    overlap = 1.0 - diff_area / tot_area if tot_area > 0 else float("nan")

    return {
        "pearson_corr": pearson,
        "cosine": cosine,
        "rmse": rmse,
        "spectral_overlap": overlap,
    }


def report_changes(changes: Dict[str, Dict]) -> str:
    """Pretty-print the highlight_changes result."""
    lines = []
    lines.append(f"{'Band':22s} {'A_before':>10s} {'A_after':>10s} "
                 f"{'Delta':>10s} {'Rel %':>8s}  Flag")
    lines.append("-" * 75)
    for name, c in changes.items():
        rel_pct = (c["rel_change"] * 100.0
                   if c["rel_change"] != float("inf") else float("inf"))
        rel_str = f"{rel_pct:+8.1f}" if rel_pct != float("inf") else "    inf"
        lines.append(
            f"{name:22s} {c['A_before']:10.4f} {c['A_after']:10.4f} "
            f"{c['delta']:+10.4f} {rel_str}  {c['flag']}"
        )
    return "\n".join(lines)
