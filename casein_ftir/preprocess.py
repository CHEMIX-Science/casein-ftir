"""
casein_ftir.preprocess
======================

Prétraitement de spectres FTIR :
    - Correction de ligne de base : ALS, polynomial, rubberband
    - Lissage : Savitzky-Golay
    - Normalisation : min-max, vector, area, par pic de référence

Toutes les fonctions retournent un NOUVEAU Spectrum (immutabilité).
"""

from __future__ import annotations
import numpy as np
from scipy.signal import savgol_filter
from .io_module import Spectrum


# ---------------------------------------------------------------------------
# Baseline correction
# ---------------------------------------------------------------------------

def baseline_als(spec: Spectrum, lam: float = 1e5, p: float = 0.001,
                 n_iter: int = 10) -> Spectrum:
    """Asymmetric Least Squares (Eilers & Boelens, 2005).

    Robust baseline for FTIR/Raman spectra.

    Parameters
    ----------
    lam : float
        Smoothness parameter. Larger = smoother baseline (1e4 - 1e7).
    p   : float
        Asymmetry parameter (0 < p < 1). 0.001 - 0.01 typical for IR.
    n_iter : int
        Number of iterations.
    """
    try:
        from pybaselines.whittaker import asls
        baseline, _ = asls(spec.absorbance, lam=lam, p=p, max_iter=n_iter)
    except ImportError:
        # Fallback : pure-numpy/scipy implementation
        baseline = _als_numpy(spec.absorbance, lam=lam, p=p, n_iter=n_iter)
    out = spec.copy()
    out.absorbance = spec.absorbance - baseline
    out.name = f"{spec.name}_baselined"
    out.metadata["baseline_method"] = f"ALS(lam={lam}, p={p})"
    return out


def _als_numpy(y: np.ndarray, lam: float, p: float, n_iter: int) -> np.ndarray:
    """Pure-numpy ALS fallback (Eilers 2005)."""
    from scipy.sparse import diags, eye, csc_matrix
    from scipy.sparse.linalg import spsolve
    L = len(y)
    if L < 3:
        raise ValueError("ALS baseline needs at least three points")
    D = diags([1.0, -2.0, 1.0], [0, -1, -2], shape=(L, L - 2))
    D = lam * D.dot(D.transpose())
    w = np.ones(L)
    for _ in range(n_iter):
        W = diags([w], [0], shape=(L, L))
        Z = csc_matrix(W + D)
        z = spsolve(Z, w * y)
        w = np.where(y > z, p, 1 - p)
    return z


def baseline_polynomial(spec: Spectrum, degree: int = 3) -> Spectrum:
    """Subtract a polynomial baseline of given degree."""
    coefs = np.polyfit(spec.wavenumber, spec.absorbance, degree)
    baseline = np.polyval(coefs, spec.wavenumber)
    out = spec.copy()
    out.absorbance = spec.absorbance - baseline
    out.name = f"{spec.name}_poly{degree}"
    out.metadata["baseline_method"] = f"polynomial(degree={degree})"
    return out


def baseline_rubberband(spec: Spectrum) -> Spectrum:
    """Rubberband baseline (lower convex hull). Simple and robust."""
    x = spec.wavenumber
    y = spec.absorbance
    # Sort ascending for the convex hull algorithm.
    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    n = len(xs)

    # Compute the lower convex hull using a monotone-chain pass.
    # We process points in ascending x; the resulting list of vertices is
    # the lower hull. This is more robust than scipy's ConvexHull, which
    # returns the *full* hull (upper + lower) without distinction.
    lower = []
    for i in range(n):
        while len(lower) >= 2:
            x1, y1 = lower[-2]
            x2, y2 = lower[-1]
            x3, y3 = xs[i], ys[i]
            # cross-product test : we want right turns only (lower hull)
            cross = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
            if cross <= 0:
                lower.pop()
            else:
                break
        lower.append((xs[i], ys[i]))

    lx = np.array([p[0] for p in lower])
    ly = np.array([p[1] for p in lower])
    baseline = np.interp(xs, lx, ly)

    # Restore original ordering
    inv = np.argsort(order)
    baseline = baseline[inv]

    out = spec.copy()
    out.absorbance = spec.absorbance - baseline
    out.name = f"{spec.name}_rubberband"
    out.metadata["baseline_method"] = "rubberband"
    return out


# ---------------------------------------------------------------------------
# Smoothing
# ---------------------------------------------------------------------------

def smooth_savgol(spec: Spectrum, window: int = 11, poly: int = 3) -> Spectrum:
    """Savitzky-Golay smoothing.

    Parameters
    ----------
    window : odd int, must be > poly. Typical: 9, 11, 15.
    poly   : polynomial order. Typical: 2 or 3.
    """
    if window % 2 == 0:
        window += 1
    if window <= poly:
        window = poly + 3 if (poly + 3) % 2 == 1 else poly + 4
    out = spec.copy()
    if window > spec.n_points:
        raise ValueError(f"Smoothing requires at least {window} points; disable smoothing for this input")
    spacing = np.diff(spec.wavenumber)
    if spacing.size and not np.allclose(spacing, spacing[0], rtol=1e-3, atol=1e-6):
        raise ValueError("Savitzky-Golay smoothing requires an evenly spaced axis; resample or disable smoothing")
    out.absorbance = savgol_filter(spec.absorbance, window, poly)
    out.name = f"{spec.name}_smoothed"
    out.metadata["smoothing"] = f"savgol(window={window}, poly={poly})"
    return out


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

def normalize_minmax(spec: Spectrum) -> Spectrum:
    """Scale absorbance to [0, 1]."""
    y = spec.absorbance
    mn, mx = float(y.min()), float(y.max())
    if mx - mn < 1e-12:
        return spec.copy()
    out = spec.copy()
    out.absorbance = (y - mn) / (mx - mn)
    out.name = f"{spec.name}_minmax"
    out.metadata["normalization"] = "min-max"
    return out


def normalize_vector(spec: Spectrum) -> Spectrum:
    """Unit L2 norm (Euclidean)."""
    y = spec.absorbance
    norm = float(np.linalg.norm(y))
    if norm < 1e-12:
        return spec.copy()
    out = spec.copy()
    out.absorbance = y / norm
    out.name = f"{spec.name}_l2norm"
    out.metadata["normalization"] = "L2 (vector)"
    return out


def normalize_area(spec: Spectrum) -> Spectrum:
    """Unit total area under curve."""
    area = float(np.trapezoid(spec.absorbance, -spec.wavenumber))  # x descending -> use -wn
    if abs(area) < 1e-12:
        return spec.copy()
    out = spec.copy()
    out.absorbance = spec.absorbance / area
    out.name = f"{spec.name}_areanorm"
    out.metadata["normalization"] = "unit area"
    return out


def normalize_to_band(spec: Spectrum, center: float, halfwidth: float = 30.0
                      ) -> Spectrum:
    """Set the max of a chosen band to 1.0.

    Useful : normalise to amide I so that ratios A_X/A_amideI become directly
    comparable across samples.
    """
    mask = (spec.wavenumber >= center - halfwidth) & \
           (spec.wavenumber <= center + halfwidth)
    if not np.any(mask):
        raise ValueError(f"No data in window around {center} cm-1")
    peak_val = float(spec.absorbance[mask].max())
    if peak_val < 1e-12:
        return spec.copy()
    out = spec.copy()
    out.absorbance = spec.absorbance / peak_val
    out.name = f"{spec.name}_norm@{center:.0f}"
    out.metadata["normalization"] = f"band {center} cm-1 -> 1.0"
    return out


# ---------------------------------------------------------------------------
# Convenience pipeline
# ---------------------------------------------------------------------------

def full_preprocess(spec: Spectrum,
                    baseline: str = "als",
                    smooth: bool = True,
                    normalize: str = "amide_I",
                    **kwargs) -> Spectrum:
    """Default preprocessing pipeline : baseline -> smooth -> normalize.

    Parameters
    ----------
    baseline : 'als', 'poly', 'rubberband', or None
    smooth   : bool, apply Savitzky-Golay
    normalize : 'minmax', 'vector', 'area', 'amide_I', or None
    """
    out = spec.copy()
    if baseline == "als":
        out = baseline_als(out)
    elif baseline == "poly":
        out = baseline_polynomial(out, degree=kwargs.get("poly_degree", 3))
    elif baseline == "rubberband":
        out = baseline_rubberband(out)
    elif baseline is None:
        pass
    else:
        raise ValueError(f"Unknown baseline method: {baseline}")

    if smooth:
        out = smooth_savgol(out, window=kwargs.get("smooth_window", 11),
                            poly=kwargs.get("smooth_poly", 3))

    if normalize == "minmax":
        out = normalize_minmax(out)
    elif normalize == "vector":
        out = normalize_vector(out)
    elif normalize == "area":
        out = normalize_area(out)
    elif normalize == "amide_I":
        out = normalize_to_band(out, center=1648.0, halfwidth=30.0)
    elif normalize is None:
        pass
    else:
        raise ValueError(f"Unknown normalization method: {normalize}")

    return out
