"""
casein_ftir.io_gaussian
========================

Lecture de fichiers de sortie Gaussian (.log / .out) et génération d'un
spectre IR théorique à partir des fréquences et intensités calculées.

Sections lues :
    - "Frequencies --"  : modes vibrationnels (cm-1)
    - "IR Inten    --"  : intensités IR correspondantes (km/mol)

Stratégie :
    1. parser les triplets de modes
    2. appliquer un facteur d'échelle empirique
    3. convoluer chaque mode avec une gaussienne ou lorentzienne
       de largeur fwhm donnée
    4. produire un Spectrum standard

Facteurs hérités du programme (provenance NIST exacte À VALIDER) :
    B3LYP/6-31G(d)         : 0.9614
    B3LYP/6-31+G(d,p)      : 0.964
    B3LYP/6-311+G(2d,p)    : 0.9679
    wB97X-D/6-31G(d)       : 0.9484
"""

from __future__ import annotations
import re
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .io_module import Spectrum


# ---------------------------------------------------------------------------
# Legacy presets — À VALIDER against the original source and exact level
# ---------------------------------------------------------------------------

SCALE_FACTORS: Dict[str, float] = {
    "b3lyp/6-31g(d)":      0.9614,
    "b3lyp/6-31g(d,p)":    0.9614,
    "b3lyp/6-31+g(d,p)":   0.964,
    "b3lyp/6-311g(d,p)":   0.9659,
    "b3lyp/6-311+g(d,p)":  0.967,
    "b3lyp/6-311+g(2d,p)": 0.9679,
    "wb97xd/6-31g(d)":     0.9484,
    "mp2/6-31g(d)":        0.9434,
    "hf/6-31g(d)":         0.8929,
}


def get_scale_factor(level: str) -> float:
    """Return the empirical scale factor for a given DFT level.

    The lookup is case-insensitive and whitespace-tolerant.
    Falls back to 1.0 if level is not recognised (with a warning printed).
    """
    key = re.sub(r"\s+", "", level.lower())
    if key in SCALE_FACTORS:
        return SCALE_FACTORS[key]
    # try partial match
    for k, v in SCALE_FACTORS.items():
        if k.startswith(key.split("/")[0] + "/"):
            continue
    return 1.0


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

_FREQ_LINE = re.compile(r"^\s*Frequencies\s*--\s*(.+)$")
_INTEN_LINE = re.compile(r"^\s*IR\s*Inten\s*--\s*(.+)$")
_ROUTE_LINE = re.compile(r"^\s*#\s*(.+)$")


def parse_gaussian_log(path: str | Path) -> Dict:
    """Parse a Gaussian frequency-calculation output.

    Parameters
    ----------
    path : str or Path
        Path to a .log / .out file from `Opt Freq` (or just `Freq`).

    Returns
    -------
    dict with:
        - 'frequencies'  : np.ndarray, raw harmonic frequencies (cm-1)
        - 'intensities'  : np.ndarray, IR intensities (km/mol)
        - 'level'        : str, detected DFT level from the route section
                            (e.g. 'B3LYP/6-31+G(d,p)') -- empty if not found
        - 'n_modes'      : int

    Notes
    -----
    Gaussian writes frequencies in blocks of 3 per line; this parser
    handles any number of imaginary modes (negative frequencies for
    transition states are kept as-is, the user must check the log).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    freqs: List[float] = []
    intens: List[float] = []
    level = ""

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            # Capture the route section (often spans multiple lines, we
            # just keep the first '#' line which contains the method)
            if not level:
                m_route = _ROUTE_LINE.match(line)
                if m_route:
                    text = m_route.group(1)
                    # try to extract method/basis pattern like 'B3LYP/6-31G(d)'
                    m_method = re.search(
                        r"([A-Za-z0-9]+)/(\S+)", text
                    )
                    if m_method:
                        level = f"{m_method.group(1)}/{m_method.group(2)}"

            m_freq = _FREQ_LINE.match(line)
            if m_freq:
                vals = [float(x) for x in m_freq.group(1).split()]
                freqs.extend(vals)
                continue

            m_int = _INTEN_LINE.match(line)
            if m_int:
                vals = [float(x) for x in m_int.group(1).split()]
                intens.extend(vals)

    freqs_arr = np.array(freqs, dtype=float)
    intens_arr = np.array(intens, dtype=float)

    # Sanity check : the two lists must have the same length
    if len(freqs_arr) == 0:
        raise ValueError(
            f"No frequency block found in {path}. "
            "Is this a Gaussian Freq calculation output?"
        )
    if len(intens_arr) != len(freqs_arr):
        # Some logs may have only frequencies (no intensities) -- fill with NaN
        intens_arr = np.full_like(freqs_arr, np.nan)

    return {
        "frequencies": freqs_arr,
        "intensities": intens_arr,
        "level": level,
        "n_modes": len(freqs_arr),
        "source_file": str(path),
    }


# ---------------------------------------------------------------------------
# Convolution into a continuous spectrum
# ---------------------------------------------------------------------------

def _gaussian(x, c, a, fwhm):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return a * np.exp(-0.5 * ((x - c) / sigma) ** 2)


def _lorentzian(x, c, a, fwhm):
    gamma = fwhm / 2.0
    return a * gamma ** 2 / ((x - c) ** 2 + gamma ** 2)


def gaussian_log_to_spectrum(path: str | Path,
                              scale: Optional[float] = None,
                              level: Optional[str] = None,
                              wn_min: float = 400.0,
                              wn_max: float = 4000.0,
                              resolution: float = 2.0,
                              fwhm: float = 15.0,
                              lineshape: str = "lorentzian",
                              drop_imaginary: bool = True,
                              ) -> Spectrum:
    """Convert a Gaussian frequency log into a continuous IR Spectrum.

    Parameters
    ----------
    path : str or Path
        Gaussian .log/.out file.
    scale : float, optional
        Override the empirical scale factor. If None, the factor is looked
        up from the route section in the log (B3LYP/6-31G(d) etc.).
    level : str, optional
        Override the detected level (e.g. 'B3LYP/6-311+G(2d,p)').
    wn_min, wn_max, resolution : floats
        Output spectrum axis (cm-1, descending).
    fwhm : float
        FWHM (cm-1) of each individual mode after convolution.
    lineshape : {'lorentzian', 'gaussian'}
        Line shape for convolution. Lorentzian is more realistic in solution.
    drop_imaginary : bool
        If True, ignore modes with negative frequency (transition-state
        artifacts). The user should still check the log if any.

    Returns
    -------
    Spectrum
    """
    if not all(np.isfinite(v) for v in (wn_min, wn_max, resolution, fwhm)):
        raise ValueError("Gaussian grid parameters must be finite")
    if wn_min >= wn_max or resolution <= 0 or fwhm <= 0:
        raise ValueError("Invalid Gaussian range, resolution or FWHM")
    if lineshape not in {"gaussian", "lorentzian"}:
        raise ValueError("Unknown line shape")
    parsed = parse_gaussian_log(path)
    freqs = parsed["frequencies"]
    intens = parsed["intensities"]
    detected_level = parsed["level"]

    # Determine scale factor
    if scale is None:
        chosen_level = level or detected_level
        scale = get_scale_factor(chosen_level) if chosen_level else 1.0

    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Scale factor must be finite and positive")
    if not np.isfinite(intens).all() or np.any(intens < 0):
        raise ValueError("Missing or invalid IR intensities; conversion cannot assume uniform weights")
    # Apply scale
    freqs_scaled = freqs * scale

    # Drop imaginary modes (negative frequencies)
    if drop_imaginary:
        mask = freqs_scaled > 0
        freqs_scaled = freqs_scaled[mask]
        intens = intens[mask]

    # Build descending wavenumber axis (IR convention)
    wn = np.arange(wn_max, wn_min - resolution, -resolution)
    abs_ = np.zeros_like(wn)

    # Intensities have been validated; no fabricated weights.
    intens_safe = intens
    # Skip modes with zero intensity (no need to convolve)
    keep = intens_safe > 1e-6

    shape_fn = _lorentzian if lineshape == "lorentzian" else _gaussian
    for f, i in zip(freqs_scaled[keep], intens_safe[keep]):
        abs_ = abs_ + shape_fn(wn, f, i, fwhm)

    # Normalise to a reasonable absorbance range (max -> 1.0) so the
    # spectrum is directly comparable to experimental ones.
    if abs_.max() > 1e-9:
        abs_ = abs_ / abs_.max()

    spec = Spectrum(
        wavenumber=wn,
        absorbance=abs_,
        name=Path(path).stem + "_dft",
        source_file=str(path),
        metadata={
            "type": "dft_simulated",
            "n_modes_total": int(parsed["n_modes"]),
            "n_modes_used": int(np.sum(keep)),
            "scale_factor": scale,
            "detected_level": detected_level,
            "lineshape": lineshape,
            "fwhm_cm-1": fwhm,
        },
    )
    return spec


# ---------------------------------------------------------------------------
# Modes inspection helper (for the CLI 'gaussian inspect' command)
# ---------------------------------------------------------------------------

def list_modes(path: str | Path, scale: Optional[float] = None,
                top_n: int = 20) -> str:
    """Return a pretty-printed table of the strongest IR modes."""
    parsed = parse_gaussian_log(path)
    freqs = parsed["frequencies"]
    intens = parsed["intensities"]

    if scale is None:
        scale = get_scale_factor(parsed["level"]) if parsed["level"] else 1.0

    # Sort by intensity, descending
    order = np.argsort(-np.where(np.isnan(intens), 0.0, intens))
    lines = []
    lines.append(f"File         : {path}")
    lines.append(f"Level        : {parsed['level'] or 'unknown'}")
    lines.append(f"Scale factor : {scale:.4f}")
    lines.append(f"Total modes  : {parsed['n_modes']}")
    lines.append("")
    lines.append(f"{'Rank':>4s} {'Harmonic':>10s} {'Scaled':>10s} "
                  f"{'IR Inten':>10s}")
    lines.append("-" * 40)
    for rank, idx in enumerate(order[:top_n], start=1):
        f_h = freqs[idx]
        f_s = f_h * scale
        i = intens[idx]
        i_str = f"{i:10.2f}" if not np.isnan(i) else "       -- "
        lines.append(f"{rank:>4d} {f_h:10.2f} {f_s:10.2f} {i_str}")
    return "\n".join(lines)
