"""
casein_ftir.report
==================

Export de rapports d'analyse en texte simple ou HTML.

Pour un rapport plus riche (PDF avec graphiques), utiliser matplotlib
en combinaison via le script CLI.
"""

from __future__ import annotations
from pathlib import Path
from datetime import datetime
from html import escape

def _e(value):
    return escape(str(value), quote=True)

CAUTION = ("PROPOSITION — Exploratory model outputs. Band assignments and indicators "
           "require independent validation (À VALIDER). Area fractions are not proven "
           "structural populations. No purity, safety, formaldehyde concentration or "
           "crosslink conversion is established by this report.")
from typing import Dict, Optional

from .io_module import Spectrum
from .database import REFERENCES


def text_report(spectrum: Spectrum,
                peaks: Optional[Dict] = None,
                ratios: Optional[Dict] = None,
                deconv: Optional[Dict] = None,
                galalithe: Optional[Dict] = None,
                changes: Optional[Dict] = None,
                reference_analysis: Optional[Dict] = None,
                output: Optional[str] = None) -> str:
    """Build a plain-text analysis report."""
    lines = []
    lines.append("=" * 75)
    lines.append("CHEMIX — CASEIN FTIR ANALYSIS REPORT")
    lines.append(CAUTION)
    lines.append("=" * 75)
    lines.append(f"Generated:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Spectrum:    {spectrum.name}")
    if spectrum.source_file:
        lines.append(f"Source file: {Path(spectrum.source_file).name}")
    wn1, wn2 = spectrum.range
    lines.append(f"Range:       {wn2:.1f} - {wn1:.1f} cm-1 ({spectrum.n_points} points)")
    if spectrum.metadata:
        lines.append("Metadata:")
        for k, v in spectrum.metadata.items():
            if k not in {"type", "baseline_method", "smoothing", "normalization", "reference_id", "scale_factor", "lineshape", "fwhm_cm-1"}:
                continue
            lines.append(f"  {k}: {v}")
    lines.append("")

    if peaks is not None:
        from .peak_detection import report_peaks
        lines.append("--- Peak Detection ---")
        lines.append(report_peaks(peaks))
        lines.append("")

    if ratios is not None:
        lines.append("--- Diagnostic Ratios ---")
        for k, v in ratios.items():
            lines.append(f"  {k:25s} = {v:.4f}")
        lines.append("")

    if deconv is not None:
        from .amide_deconv import report_deconvolution
        lines.append(report_deconvolution(deconv))
        lines.append("")

    if galalithe is not None:
        lines.append("--- Galalithe Polymerization Indicators ---")
        for k, v in galalithe.items():
            lines.append(f"  {k:30s} = {v:.4f}")
        lines.append("")

    if changes is not None:
        from .compare import report_changes
        lines.append("--- Before / After comparison ---")
        lines.append(report_changes(changes))
        lines.append("")

    if reference_analysis is not None:
        from .reference_module import report_reference_comparison
        lines.append(report_reference_comparison(reference_analysis))
        lines.append("")

    lines.append("=" * 75)
    lines.append("REFERENCES (literature used for band assignments)")
    lines.append("=" * 75)
    for key, ref in REFERENCES.items():
        lines.append(f"[{key}] {ref}")

    text = "\n".join(lines)
    if output:
        Path(output).write_text(text, encoding="utf-8")
    return text


def html_report(spectrum: Spectrum,
                peaks: Optional[Dict] = None,
                ratios: Optional[Dict] = None,
                deconv: Optional[Dict] = None,
                galalithe: Optional[Dict] = None,
                changes: Optional[Dict] = None,
                reference_analysis: Optional[Dict] = None,
                output: Optional[str] = None) -> str:
    """Build a minimal HTML report (no external dependencies)."""
    parts = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
             "<title>FTIR Analysis Report</title>",
             "<style>",
             "body{font-family:sans-serif;max-width:900px;margin:2em auto;padding:1em}",
             "table{border-collapse:collapse;margin:1em 0}",
             "td,th{padding:4px 12px;border:1px solid #999;text-align:right}",
             "th{background:#eee}",
             "h2{border-bottom:2px solid #333;margin-top:2em}",
             ".meta{color:#666;font-size:0.9em}",
             "</style></head><body>"]

    parts.append("<h1>CHEMIX — FTIR Analysis Report</h1>")
    parts.append(f"<p><strong>{_e(CAUTION)}</strong></p>")
    parts.append(f"<p class='meta'>Generated {datetime.now():%Y-%m-%d %H:%M}</p>")
    parts.append(f"<p><b>Spectrum:</b> {_e(spectrum.name)}</p>")
    if spectrum.source_file:
        parts.append(f"<p class='meta'>Source: {_e(Path(spectrum.source_file).name)}</p>")

    if peaks is not None:
        parts.append("<h2>Detected peaks</h2>")
        parts.append("<table><tr><th>Band</th><th>Position</th>"
                     "<th>Expected</th><th>Shift</th><th>Absorbance</th></tr>")
        for name, p in peaks.items():
            if p.found:
                parts.append(
                    f"<tr><td>{_e(name)}</td><td>{p.position:.1f}</td>"
                    f"<td>{p.expected_position:.1f}</td>"
                    f"<td>{p.shift:+.1f}</td>"
                    f"<td>{p.absorbance:.4f}</td></tr>"
                )
        parts.append("</table>")

    if ratios is not None:
        parts.append("<h2>Diagnostic ratios</h2><table>")
        for k, v in ratios.items():
            parts.append(f"<tr><td>{_e(k)}</td><td>{v:.4f}</td></tr>")
        parts.append("</table>")

    if deconv is not None:
        parts.append("<h2>Amide I deconvolution</h2>")
        parts.append(f"<p>Fit R² = {deconv['fit_quality_r2']:.4f}, "
                     f"total area = {deconv['total_area']:.4f}</p>")
        parts.append("<table><tr><th>Sub-band</th><th>Center</th>"
                     "<th>FWHM</th><th>Area</th><th>%</th></tr>")
        for name, comp in deconv["components"].items():
            pct = deconv["structure_pct_per_subband"][name]
            parts.append(
                f"<tr><td>{_e(name)}</td><td>{comp['center']:.1f}</td>"
                f"<td>{comp['fwhm']:.2f}</td>"
                f"<td>{comp['area']:.4f}</td><td>{pct:.2f}</td></tr>"
            )
        parts.append("</table>")
        parts.append("<h3>Model-assigned area fractions</h3><table>")
        for k, v in deconv["summary"].items():
            parts.append(f"<tr><td>{_e(k)}</td><td>{v:.2f}%</td></tr>")
        parts.append("</table>")

    if galalithe is not None:
        parts.append("<h2>Galalithe indicators</h2><table>")
        for k, v in galalithe.items():
            parts.append(f"<tr><td>{_e(k)}</td><td>{v:.4f}</td></tr>")
        parts.append("</table>")

    if changes is not None:
        parts.append("<h2>Before / After comparison</h2>")
        parts.append("<table><tr><th>Band</th><th>A before</th><th>A after</th>"
                     "<th>Δ</th><th>Rel %</th><th>Flag</th></tr>")
        for name, c in changes.items():
            rel = c["rel_change"]
            rel_str = f"{rel*100:+.1f}" if rel != float("inf") else "∞"
            parts.append(
                f"<tr><td>{_e(name)}</td><td>{c['A_before']:.4f}</td>"
                f"<td>{c['A_after']:.4f}</td>"
                f"<td>{c['delta']:+.4f}</td>"
                f"<td>{rel_str}</td><td><b>{_e(c['flag'])}</b></td></tr>"
            )
        parts.append("</table>")

    if reference_analysis is not None:
        sim = reference_analysis["similarity"]
        parts.append("<h2>Reference comparison</h2>")
        parts.append(f"<p><b>Heuristic similarity category :</b> "
                       f"{_e(reference_analysis['overall_quality'].upper())}</p>")
        parts.append(f"<p><i>{_e(reference_analysis['recommended_action'])}</i></p>")
        parts.append("<table>")
        parts.append(f"<tr><td>Pearson correlation</td>"
                       f"<td>{sim['pearson_corr']:.4f}</td></tr>")
        parts.append(f"<tr><td>Cosine similarity</td>"
                       f"<td>{sim['cosine']:.4f}</td></tr>")
        parts.append(f"<tr><td>RMSE</td><td>{sim['rmse']:.4f}</td></tr>")
        parts.append(f"<tr><td>Spectral overlap</td>"
                       f"<td>{sim['spectral_overlap']:.4f}</td></tr>")
        parts.append("</table>")
        parts.append("<h3>Per-band deviation</h3>")
        parts.append("<table><tr><th>Band</th><th>Sample area</th>"
                       "<th>Ref area</th><th>Diff %</th>"
                       "<th>Shift cm⁻¹</th><th>Flag</th></tr>")
        for name, d in reference_analysis["band_deviations"].items():
            diff = d["rel_area_diff_pct"]
            shift = d["shift_cm-1"]
            import math
            diff_s = (f"{diff:+.1f}" if diff not in (float("inf"),)
                       and not math.isnan(diff) else "—")
            shift_s = f"{shift:+.1f}" if not math.isnan(shift) else "—"
            parts.append(
                f"<tr><td>{_e(name)}</td><td>{d['sample_area']:.4f}</td>"
                f"<td>{d['reference_area']:.4f}</td>"
                f"<td>{diff_s}</td><td>{shift_s}</td>"
                f"<td><b>{_e(d['flag'])}</b></td></tr>"
            )
        parts.append("</table>")
        # Peak-by-peak matching table (if available)
        if "peak_matching" in reference_analysis:
            pm = reference_analysis["peak_matching"]
            parts.append("<h3>Peak-by-peak matching</h3>")
            parts.append(
                f"<p><b>Matched : {pm['n_matched']}/{pm['n_ref_peaks']}</b> "
                f"reference peaks (within ±{pm['tolerance_cm-1']:.0f} cm⁻¹).<br>"
                f"Match rate : {100*pm['match_rate']:.1f}% &nbsp; "
                f"Intensity-weighted : "
                f"{100*pm['weighted_match_rate']:.1f}%</p>"
            )
            parts.append("<table><tr><th>Ref peak (cm⁻¹)</th><th>Ref I</th>"
                           "<th>Nearest sample (cm⁻¹)</th><th>Sample I</th>"
                           "<th>Δ (cm⁻¹)</th><th>Match</th></tr>")
            for m in pm["matches"]:
                import math
                if math.isnan(m["delta_wn"]):
                    parts.append(
                        f"<tr><td>{m['ref_wn']:.1f}</td>"
                        f"<td>{m['ref_intensity']:.4f}</td>"
                        f"<td>—</td><td>—</td><td>—</td>"
                        f"<td><b>no</b></td></tr>"
                    )
                else:
                    flag = "<b>YES</b>" if m["matched"] else "no"
                    parts.append(
                        f"<tr><td>{m['ref_wn']:.1f}</td>"
                        f"<td>{m['ref_intensity']:.4f}</td>"
                        f"<td>{m['nearest_sample_wn']:.1f}</td>"
                        f"<td>{m['nearest_sample_intensity']:.4f}</td>"
                        f"<td>{m['delta_wn']:+.1f}</td>"
                        f"<td>{flag}</td></tr>"
                    )
            parts.append("</table>")

    parts.append("<h2>References</h2><ol>")
    for ref in REFERENCES.values():
        parts.append(f"<li>{_e(ref)}</li>")
    parts.append("</ol>")
    parts.append("</body></html>")

    html = "\n".join(parts)
    if output:
        Path(output).write_text(html, encoding="utf-8")
    return html
