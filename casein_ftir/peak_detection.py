"""
casein_ftir.peak_detection
==========================

Détection automatique des pics caractéristiques dans des fenêtres
prédéfinies (cf. casein_ftir.database).

Pour chaque bande attendue, on cherche le maximum local dans sa fenêtre
de tolérance, en utilisant scipy.signal.find_peaks pour gérer les
épaulements et les pics multiples.
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional
from scipy.signal import find_peaks

from .io_module import Spectrum
from .database import CASEIN_BANDS, GALALITHE_MARKERS, IRBand


@dataclass
class DetectedPeak:
    """Result of one peak detection."""
    name: str                    # band name (e.g. 'amide_I')
    position: float              # cm^-1
    absorbance: float            # local max value
    expected_position: float     # theoretical center (from database)
    shift: float                 # position - expected_position (cm^-1)
    window: tuple                # (min, max) cm^-1
    found: bool                  # True if a peak was located

    def as_dict(self) -> dict:
        return asdict(self)


def _local_max_in_window(spec: Spectrum, wn_min: float, wn_max: float,
                          prominence: Optional[float] = None) -> Optional[tuple]:
    """Return (wavenumber, absorbance) of the largest peak in window."""
    mask = (spec.wavenumber >= wn_min) & (spec.wavenumber <= wn_max)
    if not np.any(mask):
        return None
    x = spec.wavenumber[mask]
    y = spec.absorbance[mask]

    # IR convention : wn descending. Ensure ascending for find_peaks.
    if x[0] > x[-1]:
        x = x[::-1]
        y = y[::-1]

    if prominence is None:
        prominence = 0.01 * float(np.ptp(y))

    peaks, props = find_peaks(y, prominence=prominence)
    if len(peaks) > 0:
        # Pick the most prominent peak
        idx = peaks[np.argmax(props["prominences"])]
        return float(x[idx]), float(y[idx])

    # No real peak found -> return the absolute max in the window
    # (useful for broad bands like amide A)
    idx = int(np.argmax(y))
    return float(x[idx]), float(y[idx])


def detect_bands(spec: Spectrum,
                 band_set: str = "casein",
                 prominence: Optional[float] = None) -> Dict[str, DetectedPeak]:
    """Detect all expected bands in a spectrum.

    Parameters
    ----------
    spec : Spectrum
    band_set : {'casein', 'galalithe', 'all'}
    prominence : float, optional
        Minimum prominence in absorbance units. Auto-scaled if None.
    """
    if band_set == "casein":
        bands = CASEIN_BANDS
    elif band_set == "galalithe":
        bands = GALALITHE_MARKERS
    elif band_set == "all":
        bands = {**CASEIN_BANDS, **GALALITHE_MARKERS}
    else:
        raise ValueError(f"Unknown band_set: {band_set!r}")

    results: Dict[str, DetectedPeak] = {}
    for name, band in bands.items():
        wn_min, wn_max = band.window
        res = _local_max_in_window(spec, wn_min, wn_max, prominence=prominence)
        if res is None:
            results[name] = DetectedPeak(
                name=name, position=np.nan, absorbance=np.nan,
                expected_position=band.center, shift=np.nan,
                window=band.window, found=False,
            )
        else:
            pos, ab = res
            results[name] = DetectedPeak(
                name=name, position=pos, absorbance=ab,
                expected_position=band.center,
                shift=pos - band.center,
                window=band.window, found=True,
            )
    return results


def report_peaks(peaks: Dict[str, DetectedPeak]) -> str:
    """Pretty-print the peak detection results."""
    lines = []
    lines.append(f"{'Band':20s} {'Position':>10s} {'Expected':>10s} "
                 f"{'Shift':>8s} {'Abs':>10s} {'Found':>6s}")
    lines.append("-" * 70)
    for name, p in peaks.items():
        if p.found:
            lines.append(
                f"{name:20s} {p.position:10.1f} {p.expected_position:10.1f} "
                f"{p.shift:+8.1f} {p.absorbance:10.4f} {'YES':>6s}"
            )
        else:
            lines.append(
                f"{name:20s} {'--':>10s} {p.expected_position:10.1f} "
                f"{'--':>8s} {'--':>10s} {'NO':>6s}"
            )
    return "\n".join(lines)
