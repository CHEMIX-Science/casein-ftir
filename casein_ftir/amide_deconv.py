"""
casein_ftir.amide_deconv
========================

Déconvolution gaussienne de la bande amide I (1600-1700 cm-1) pour
l'analyse de la structure secondaire.

Six sous-bandes selon Barth (2007) et Sadat & Joye (2020) :
    - beta_aggregate   (~1618 cm-1)
    - beta_sheet       (~1632 cm-1)
    - random_coil      (~1644 cm-1)
    - alpha_helix      (~1654 cm-1)
    - beta_turn        (~1670 cm-1)
    - beta_antiparallel(~1685 cm-1)

PROPOSITION — Les fractions sont des aires relatives de composantes du modèle.
À VALIDER — Leur conversion en populations structurales nécessite une validation externe.
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple, Optional
from .io_module import Spectrum
from .database import AMIDE_I_SUBBANDS


def _gaussian(x, c, a, fwhm):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return a * np.exp(-0.5 * ((x - c) / sigma) ** 2)


def deconvolve_amide_I(spec: Spectrum,
                       wn_min: float = 1600.0,
                       wn_max: float = 1700.0,
                       baseline_correct: bool = True,
                       default_fwhm: float = 18.0,
                       max_fwhm: float = 35.0,
                       min_fwhm: float = 8.0) -> Dict:
    """Decompose the amide I band into Gaussian sub-bands.

    Uses lmfit if available, falls back to scipy.optimize.curve_fit
    otherwise.

    Returns
    -------
    dict with:
        - components       : dict {name: {'center', 'amplitude', 'fwhm', 'area'}}
        - structure_pct    : dict {name: percentage of total area}
        - total_area       : float
        - fit_quality      : R^2 of fit
        - residual         : np.ndarray (residual after subtraction)
        - x, y_fit         : np.ndarrays for plotting
    """
    try:
        from lmfit.models import GaussianModel
        from lmfit import Parameters
        _has_lmfit = True
    except ImportError:
        _has_lmfit = False

    sub = spec.slice(wn_min, wn_max)
    x = sub.wavenumber.copy()
    y = sub.absorbance.copy()
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if len(x) < 19:
        raise ValueError("Amide I fitting needs at least 19 points in the selected range")
    if np.ptp(y) <= 1e-12:
        raise ValueError("Amide I fitting needs a nonconstant signal")
    # Local linear baseline subtraction
    if baseline_correct and len(x) > 1:
        y_base = np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
        y = y - y_base
        y = np.clip(y, 0.0, None)

    sub_band_names = list(AMIDE_I_SUBBANDS.keys())

    if _has_lmfit:
        return _deconvolve_lmfit(x, y, sub_band_names, default_fwhm,
                                  min_fwhm, max_fwhm)
    else:
        return _deconvolve_scipy(x, y, sub_band_names, default_fwhm,
                                  min_fwhm, max_fwhm)


def _deconvolve_lmfit(x, y, sub_band_names, default_fwhm, min_fwhm, max_fwhm):
    from lmfit.models import GaussianModel
    composite = None
    params = None
    for name in sub_band_names:
        band = AMIDE_I_SUBBANDS[name]
        prefix = f"{name}_"
        g = GaussianModel(prefix=prefix)
        if composite is None:
            composite = g
            params = g.make_params()
        else:
            composite = composite + g
            params.update(g.make_params())
        wmin, wmax = band.window
        params[f"{prefix}center"].set(value=band.center, min=wmin, max=wmax)
        sig_init = default_fwhm / 2.3548
        sig_min  = min_fwhm     / 2.3548
        sig_max  = max_fwhm     / 2.3548
        params[f"{prefix}sigma"].set(value=sig_init, min=sig_min, max=sig_max)
        params[f"{prefix}amplitude"].set(
            value=max(float(y.max()), 1e-4) * (band.rel_intensity_ref or 0.1)
                  * sig_init * np.sqrt(2 * np.pi),
            min=0.0,
        )
    result = composite.fit(y, params, x=x, method="least_squares", max_nfev=5000)

    comp_dict = {}
    total_area = 0.0
    for name in sub_band_names:
        prefix = f"{name}_"
        c   = float(result.params[f"{prefix}center"].value)
        sig = float(result.params[f"{prefix}sigma"].value)
        amp = float(result.params[f"{prefix}amplitude"].value)  # already area
        fwhm = 2.3548 * sig
        height = amp / (sig * np.sqrt(2 * np.pi))
        comp_dict[name] = {"center": c, "amplitude_height": height,
                            "fwhm": fwhm, "area": amp}
        total_area += amp

    if not result.success:
        raise RuntimeError(f"Deconvolution did not converge: {result.message}")
    return _finalize_deconv(x, y, result.best_fit, comp_dict, total_area)


def _deconvolve_scipy(x, y, sub_band_names, default_fwhm, min_fwhm, max_fwhm):
    """Pure-scipy fallback when lmfit is not installed."""
    from scipy.optimize import curve_fit

    n = len(sub_band_names)

    # Build initial guesses, bounds
    p0 = []
    lb = []
    ub = []
    for name in sub_band_names:
        band = AMIDE_I_SUBBANDS[name]
        wmin, wmax = band.window
        p0.extend([band.center,
                    max(float(y.max()) * 0.1 *
                         (band.rel_intensity_ref or 0.1), 1e-4),
                    default_fwhm])
        lb.extend([wmin, 0.0, min_fwhm])
        ub.extend([wmax, np.inf, max_fwhm])

    def model(x_, *params):
        out = np.zeros_like(x_)
        for i in range(n):
            c, amp, fwhm = params[3*i: 3*i+3]
            out = out + _gaussian(x_, c, amp, fwhm)
        return out

    try:
        popt, _ = curve_fit(model, x, y, p0=p0, bounds=(lb, ub),
                              maxfev=20000)
    except Exception as e:
        raise RuntimeError(f"Deconvolution fit failed: {e}")

    comp_dict = {}
    total_area = 0.0
    for i, name in enumerate(sub_band_names):
        c, amp, fwhm = popt[3*i: 3*i+3]
        sigma = fwhm / 2.3548
        area = amp * sigma * np.sqrt(2 * np.pi)  # area of Gaussian
        comp_dict[name] = {"center": float(c),
                            "amplitude_height": float(amp),
                            "fwhm": float(fwhm),
                            "area": float(area)}
        total_area += area

    y_fit = model(x, *popt)
    return _finalize_deconv(x, y, y_fit, comp_dict, total_area)


def _finalize_deconv(x, y, y_fit, comp_dict, total_area):
    structure_pct = {
        name: 100.0 * comp["area"] / total_area if total_area > 0 else float("nan")
        for name, comp in comp_dict.items()
    }
    summary = {
        "alpha_helix_pct": structure_pct.get("alpha_helix", 0.0),
        "beta_sheet_pct": (
            structure_pct.get("beta_sheet", 0.0)
            + structure_pct.get("beta_antiparallel", 0.0)
            + structure_pct.get("beta_aggregate", 0.0)
        ),
        "random_coil_pct": structure_pct.get("random_coil", 0.0),
        "beta_turn_pct":   structure_pct.get("beta_turn", 0.0),
    }
    ss_res = float(np.sum((y - y_fit) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "components": comp_dict,
        "structure_pct_per_subband": structure_pct,
        "summary": summary,
        "total_area": total_area,
        "fit_quality_r2": r2,
        "residual": y - y_fit,
        "x": x,
        "y_corrected": y,
        "y_fit": y_fit,
    }


def report_deconvolution(result: Dict) -> str:
    """Pretty-print the deconvolution result."""
    lines = []
    lines.append("=== AMIDE I DECONVOLUTION ===")
    lines.append(f"Fit quality R^2 = {result['fit_quality_r2']:.4f}")
    lines.append(f"Total amide I area = {result['total_area']:.4f}")
    lines.append("")
    lines.append(f"{'Sub-band':20s} {'Center':>8s} {'FWHM':>8s} {'Area':>10s} {'%':>7s}")
    lines.append("-" * 60)
    for name, comp in result["components"].items():
        pct = result["structure_pct_per_subband"][name]
        lines.append(
            f"{name:20s} {comp['center']:8.1f} {comp['fwhm']:8.2f} "
            f"{comp['area']:10.4f} {pct:6.2f}%"
        )
    lines.append("")
    lines.append("Model-assigned area fractions (not validated structural populations):")
    summary = result["summary"]
    for k, v in summary.items():
        lines.append(f"  {k:25s} {v:6.2f}%")
    return "\n".join(lines)


def export_deconvolution_csv(result: Dict, path: str) -> None:
    """Export the amide I deconvolution to a CSV file ready for plotting.

    Columns
    -------
        wavenumber          : x-axis (cm-1), ascending
        observed            : baseline-corrected experimental absorbance
        fit_total           : sum of all 6 fitted Gaussians (the model)
        beta_aggregate      : individual sub-band at ~1618 cm-1
        beta_sheet          : individual sub-band at ~1632 cm-1
        random_coil         : individual sub-band at ~1644 cm-1
        alpha_helix         : individual sub-band at ~1654 cm-1
        beta_turn           : individual sub-band at ~1670 cm-1
        beta_antiparallel   : individual sub-band at ~1685 cm-1
        residual            : observed - fit_total

    Use this file to plot the canonical amide I deconvolution figure
    (observed envelope + sub-bands + total fit) in Excel, Origin,
    matplotlib, etc.
    """
    from pathlib import Path
    x = result["x"]
    y_obs = result["y_corrected"]
    y_fit = result["y_fit"]
    components = result["components"]

    # Compute each sub-band individually on the same x-axis
    sub_band_order = list(components.keys())
    sub_band_curves = {}
    for name in sub_band_order:
        comp = components[name]
        height = comp.get("amplitude_height")
        if height is None:
            # Recover height from area when 'amplitude_height' is missing
            sigma = comp["fwhm"] / 2.3548
            height = comp["area"] / (sigma * np.sqrt(2 * np.pi))
        sub_band_curves[name] = _gaussian(
            x, comp["center"], height, comp["fwhm"]
        )

    residual = y_obs - y_fit

    # Write CSV
    p = Path(path)
    with p.open("w", encoding="utf-8") as f:
        f.write(f"# Amide I deconvolution export\n")
        f.write(f"# Fit quality R^2 = {result['fit_quality_r2']:.6f}\n")
        f.write(f"# Total area = {result['total_area']:.6f}\n")
        f.write(f"# Sub-bands ({len(sub_band_order)}):\n")
        for name in sub_band_order:
            c = components[name]
            pct = result["structure_pct_per_subband"][name]
            f.write(f"#   {name:20s} center={c['center']:7.2f} "
                     f"FWHM={c['fwhm']:5.2f} area={c['area']:.4f} "
                     f"pct={pct:.2f}%\n")
        # Header
        header = (["wavenumber", "observed", "fit_total"]
                   + list(sub_band_order) + ["residual"])
        f.write(",".join(header) + "\n")
        # Data
        for i in range(len(x)):
            row = [f"{x[i]:.4f}", f"{y_obs[i]:.6f}", f"{y_fit[i]:.6f}"]
            for name in sub_band_order:
                row.append(f"{sub_band_curves[name][i]:.6f}")
            row.append(f"{residual[i]:.6f}")
            f.write(",".join(row) + "\n")
