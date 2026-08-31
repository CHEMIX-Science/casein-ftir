"""
casein_ftir.reference_module
=============================

Gestion d'un spectre de référence externe pour comparaison exploratoire.

Cas d'usage typique :
    - vous mesurez la caséine extraite de lait périmé (échantillon)
    - vous chargez un spectre de référence (caséine commerciale Sigma
      ou spectre théorique simulé Gaussian)
    - le module compare les deux : similarité globale, déviations par
      bande, soustraction numérique optionnelle (pas une décontamination)
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Optional, List

from .io_module import Spectrum, load_spectrum
from .preprocess import full_preprocess
from .compare import (resample_to_common_axis, similarity_metric,
                       highlight_changes)
from .quantification import integrate_band
from .database import CASEIN_BANDS, GALALITHE_MARKERS


# ---------------------------------------------------------------------------
# Reference loader (with same preprocessing as sample if requested)
# ---------------------------------------------------------------------------

def load_reference(path: str,
                    preprocess: bool = True,
                    baseline: str = "als",
                    smooth: bool = True,
                    normalize: str = "amide_I") -> Spectrum:
    """Load an external reference spectrum (DFT-simulated, commercial,
    literature, ...).

    Parameters
    ----------
    path : str
        Path to the reference spectrum (any supported format).
    preprocess : bool
        If True, apply baseline+smooth+normalize to the reference. Set
        False if the reference is already clean (e.g., theoretical).
    baseline, smooth, normalize : forwarded to full_preprocess.

    Returns
    -------
    Spectrum (with metadata flag 'is_reference' = True)
    """
    ref = load_spectrum(path)
    ref.metadata["is_reference"] = True
    if preprocess:
        ref = full_preprocess(ref, baseline=baseline, smooth=smooth,
                                  normalize=normalize)
        ref.metadata["is_reference"] = True
        ref.name = f"{ref.name}_ref"
    return ref


# ---------------------------------------------------------------------------
# Comparison sample vs reference
# ---------------------------------------------------------------------------

def compare_to_reference(sample: Spectrum, reference: Spectrum,
                            band_change_threshold: float = 0.10
                            ) -> Dict:
    """Compare a sample spectrum to a reference spectrum.

    Returns
    -------
    dict with:
        - similarity         : dict of global similarity metrics
        - band_deviations    : dict {band_name: {sample_pos, ref_pos,
                               sample_area, ref_area, abs_shift_cm-1,
                               rel_area_diff_pct, flag}}
        - overall_quality    : str ('excellent' / 'good' / 'fair' / 'poor')
                               based on Pearson correlation
        - recommended_action : human-readable text
    """
    sim = similarity_metric(sample, reference)
    pearson = sim["pearson_corr"]
    cosine = sim["cosine"]

    # Per-band deviations
    band_dev = {}
    all_bands = {**CASEIN_BANDS, **GALALITHE_MARKERS}
    for name, band in all_bands.items():
        a_sample = integrate_band(sample, *band.window)
        a_ref    = integrate_band(reference, *band.window)
        if a_ref > 1e-9:
            rel_diff = 100.0 * (a_sample - a_ref) / a_ref
        else:
            rel_diff = float("nan") if a_sample < 1e-9 else float("inf")

        # Peak positions for shift
        from .peak_detection import detect_bands
        peaks_s = detect_bands(sample.slice(band.window[0]-50,
                                              band.window[1]+50),
                                 band_set="all")
        peaks_r = detect_bands(reference.slice(band.window[0]-50,
                                                  band.window[1]+50),
                                 band_set="all")
        shift = float("nan")
        if name in peaks_s and name in peaks_r:
            if peaks_s[name].found and peaks_r[name].found:
                shift = peaks_s[name].position - peaks_r[name].position

        # Flag
        if np.isnan(rel_diff):
            flag = "n/a"
        elif rel_diff == float("inf"):
            flag = "sample-only"
        elif abs(rel_diff) < 100.0 * band_change_threshold:
            flag = "match"
        elif rel_diff > 0:
            flag = "excess"
        else:
            flag = "deficit"

        band_dev[name] = {
            "sample_area":      a_sample,
            "reference_area":   a_ref,
            "rel_area_diff_pct": rel_diff,
            "shift_cm-1":       shift,
            "flag":             flag,
        }

    # Overall quality scoring
    # If the reference is a DFT fragment (small molecule), low global
    # correlation is expected -- the fragment captures local crosslink
    # signatures, not the full protein backbone. Adapt the message.
    is_dft_fragment = (
        reference.metadata.get("type") == "dft_packaged_reference"
        or reference.metadata.get("type") == "dft_simulated"
    )

    if is_dft_fragment:
        # Don't grade overall quality -- focus on individual band markers
        # that the DFT model predicts.
        quality = "fragment_reference"
        # Identify which DFT-strong bands are detected in the sample
        # by looking at where the reference itself is strong.
        action = (
            "DFT fragment comparison only. Band coincidences do not establish "
            "crosslink identity or material composition. Inspect assignments, "
            "preprocessing and independent evidence; global similarity is not a validation."
        )
    elif pearson > 0.95:
        quality = "excellent"
        action = "High spectral similarity only; no conclusion about purity, safety or suitability."
    elif pearson > 0.85:
        quality = "good"
        action = ("Sample close to reference, minor deviations. Check "
                  "individual band flags and acquisition conditions; differences are not compound identifications.")
    elif pearson > 0.70:
        quality = "fair"
        action = ("Moderate spectral similarity. Check preprocessing, acquisition conditions "
                  "and reference relevance before interpreting differences.")
    else:
        quality = "poor"
        action = ("Sample is very different from the reference. Check "
                  "for extraction errors or wrong reference assignment.")

    return {
        "similarity":      sim,
        "band_deviations": band_dev,
        "overall_quality": quality,
        "recommended_action": action,
        "peak_matching":   peak_matching_score(sample, reference,
                                                  tolerance=15.0),
    }


def subtract_reference(sample: Spectrum, reference: Spectrum,
                         scale: float = 1.0) -> Spectrum:
    """Subtract a (scaled) reference from the sample.

    Useful to highlight what's *different* from the reference :
    contamination, additional bands, modified secondary structure.

    Parameters
    ----------
    scale : float
        Multiplicative factor applied to the reference before subtraction.
        Default 1.0 assumes both spectra are normalised (e.g., on amide I).
    """
    s, r = resample_to_common_axis(sample, reference)
    out = Spectrum(
        wavenumber=s.wavenumber.copy(),
        absorbance=s.absorbance - scale * r.absorbance,
        name=f"({s.name})-{scale:.2f}*({r.name})",
        metadata={"type": "subtracted",
                   "scale_factor": scale,
                   "reference_name": r.name},
    )
    return out


# ---------------------------------------------------------------------------
# Peak-by-peak matching (better metric for fragment-vs-material comparisons)
# ---------------------------------------------------------------------------

def peak_matching_score(sample: Spectrum, reference: Spectrum,
                          tolerance: float = 15.0,
                          ref_prominence: float = 0.05,
                          ref_height: float = 0.1,
                          sample_prominence: float = 0.05,
                          sample_height: float = 0.1,
                          ) -> Dict:
    """Peak-by-peak matching score between a sample and a reference.

    Plus pertinent que Pearson lorsque le reference est un fragment DFT
    (chimie locale) et que le sample est un matériau ou protéine complète :
    Pearson moyennerait sur des régions où l'un est plat et l'autre actif,
    ce qui n'a pas de sens. Ici on fait du *peak matching* :
        - On extrait les pics du reference avec leurs positions/intensités
        - Pour chaque pic du reference, on cherche le pic le plus proche
          dans le sample
        - On le marque "matched" si la distance < tolerance (cm-1)
        - Le score global est le % de pics du reference appariés

    Parameters
    ----------
    sample, reference : Spectrum
    tolerance : float
        Maximum cm-1 deviation for a match (default 15).
    ref_prominence, ref_height : float
        Peak-detection thresholds for the reference spectrum
        (typically a DFT reference; defaults work well after the
        max-normalisation our parser applies).
    sample_prominence, sample_height : float
        Peak-detection thresholds for the sample (experimental spectrum).

    Returns
    -------
    dict with:
        - n_ref_peaks         : int, peaks found in reference
        - n_sample_peaks      : int, peaks found in sample
        - n_matched           : int
        - match_rate          : float in [0,1]
        - matches             : list of dicts (one per reference peak):
              {ref_wn, ref_intensity, nearest_sample_wn,
               nearest_sample_intensity, delta_wn, matched (bool)}
        - tolerance_cm-1      : float (the threshold used)
        - weighted_match_rate : float, intensity-weighted match rate
    """
    from scipy.signal import find_peaks

    a, b = resample_to_common_axis(sample, reference)
    wn = a.wavenumber

    # ---- Detect peaks in both spectra ------------------------------------
    ref_idx, _ = find_peaks(b.absorbance,
                              prominence=ref_prominence,
                              height=ref_height)
    sam_idx, _ = find_peaks(a.absorbance,
                              prominence=sample_prominence,
                              height=sample_height)

    ref_wn = wn[ref_idx]
    ref_int = b.absorbance[ref_idx]
    sam_wn = wn[sam_idx]
    sam_int = a.absorbance[sam_idx]

    # ---- Match each reference peak to nearest sample peak ----------------
    matches: List[Dict] = []
    n_matched = 0
    weighted_total = 0.0
    weighted_matched = 0.0

    for w, i in zip(ref_wn, ref_int):
        if len(sam_wn) == 0:
            matches.append({
                "ref_wn": float(w), "ref_intensity": float(i),
                "nearest_sample_wn": float("nan"),
                "nearest_sample_intensity": float("nan"),
                "delta_wn": float("nan"),
                "matched": False,
            })
            weighted_total += i
            continue
        distances = np.abs(sam_wn - w)
        j = int(np.argmin(distances))
        d = float(distances[j])
        matched = d <= tolerance
        matches.append({
            "ref_wn": float(w),
            "ref_intensity": float(i),
            "nearest_sample_wn": float(sam_wn[j]),
            "nearest_sample_intensity": float(sam_int[j]),
            "delta_wn": d,
            "matched": matched,
        })
        if matched:
            n_matched += 1
            weighted_matched += i
        weighted_total += i

    match_rate = n_matched / len(ref_wn) if len(ref_wn) > 0 else float("nan")
    weighted_match_rate = (weighted_matched / weighted_total
                            if weighted_total > 0 else float("nan"))

    return {
        "n_ref_peaks":         int(len(ref_wn)),
        "n_sample_peaks":      int(len(sam_wn)),
        "n_matched":           int(n_matched),
        "match_rate":          float(match_rate),
        "weighted_match_rate": float(weighted_match_rate),
        "matches":             matches,
        "tolerance_cm-1":      float(tolerance),
    }


def report_peak_matching(result: Dict) -> str:
    """Pretty-print peak_matching_score output."""
    lines = []
    lines.append("=" * 75)
    lines.append("PEAK-BY-PEAK MATCHING (fragment DFT vs experimental)")
    lines.append("=" * 75)
    lines.append("")
    lines.append(f"Reference peaks found  : {result['n_ref_peaks']}")
    lines.append(f"Sample peaks found     : {result['n_sample_peaks']}")
    lines.append(f"Matched (within ±{result['tolerance_cm-1']:.0f} cm-1) : "
                  f"{result['n_matched']} / {result['n_ref_peaks']}")
    lines.append(f"Match rate             : {100*result['match_rate']:.1f}%")
    lines.append(f"Intensity-weighted rate: "
                  f"{100*result['weighted_match_rate']:.1f}%")
    lines.append("")
    lines.append(f"{'Ref peak':>10s} {'Ref I':>9s} | "
                  f"{'Nearest':>10s} {'Sample I':>9s} {'Delta':>8s}  Match")
    lines.append("-" * 70)
    for m in result["matches"]:
        if np.isnan(m["delta_wn"]):
            row = (f"{m['ref_wn']:10.1f} {m['ref_intensity']:9.4f} | "
                    f"{'  --  ':>10s} {'  --  ':>9s} {'  -- ':>8s}     no")
        else:
            flag = " YES " if m["matched"] else " no  "
            row = (f"{m['ref_wn']:10.1f} {m['ref_intensity']:9.4f} | "
                    f"{m['nearest_sample_wn']:10.1f} "
                    f"{m['nearest_sample_intensity']:9.4f} "
                    f"{m['delta_wn']:+8.1f}    {flag}")
        lines.append(row)
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Pretty-print global comparison
# ---------------------------------------------------------------------------

def report_reference_comparison(result: Dict) -> str:
    """Format the reference comparison result as plain text."""
    lines = []
    lines.append("=" * 75)
    lines.append("REFERENCE COMPARISON")
    lines.append("=" * 75)
    sim = result["similarity"]
    lines.append("")
    lines.append("Global similarity:")
    lines.append(f"  Pearson correlation : {sim['pearson_corr']:.4f}")
    lines.append(f"  Cosine similarity   : {sim['cosine']:.4f}")
    lines.append(f"  RMSE                : {sim['rmse']:.4f}")
    lines.append(f"  Spectral overlap    : {sim['spectral_overlap']:.4f}")
    lines.append("")
    lines.append(f"Heuristic similarity category : {result['overall_quality'].upper()}")
    lines.append(f"  -> {result['recommended_action']}")
    lines.append("")
    lines.append("Per-band deviations (sample vs reference):")
    lines.append(f"{'Band':22s} {'Sample':>10s} {'Ref':>10s} "
                  f"{'Diff %':>10s} {'Shift':>8s}  Flag")
    lines.append("-" * 75)
    for name, d in result["band_deviations"].items():
        diff = d["rel_area_diff_pct"]
        shift = d["shift_cm-1"]
        diff_str = (f"{diff:+10.1f}" if diff not in (float("inf"), float("nan"))
                       and not np.isnan(diff) else "       --")
        shift_str = (f"{shift:+8.1f}" if not np.isnan(shift) else "      --")
        lines.append(
            f"{name:22s} {d['sample_area']:10.4f} {d['reference_area']:10.4f} "
            f"{diff_str} {shift_str}  {d['flag']}"
        )
    lines.append("")
    # Append peak matching subsection if present
    if "peak_matching" in result:
        lines.append(report_peak_matching(result["peak_matching"]))
    return "\n".join(lines)
