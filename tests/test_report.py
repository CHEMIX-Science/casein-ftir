"""
Tests for casein_ftir.report
"""

import re
import pytest

from casein_ftir.peak_detection import detect_bands
from casein_ftir.quantification import diagnostic_ratios, galalithe_indicators
from casein_ftir.amide_deconv import deconvolve_amide_I
from casein_ftir.compare import highlight_changes
from casein_ftir.report import text_report, html_report


class TestTextReport:

    def test_minimal_report(self, clean_casein_spectrum):
        text = text_report(clean_casein_spectrum)
        assert isinstance(text, str)
        assert clean_casein_spectrum.name in text
        assert "REFERENCES" in text

    def test_report_with_peaks(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="casein")
        text = text_report(clean_casein_spectrum, peaks=peaks)
        assert "Peak Detection" in text
        assert "amide_I" in text

    def test_full_report_sections(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="all")
        ratios = diagnostic_ratios(clean_casein_spectrum)
        gal = galalithe_indicators(clean_casein_spectrum)
        deconv = deconvolve_amide_I(clean_casein_spectrum)
        text = text_report(clean_casein_spectrum,
                              peaks=peaks, ratios=ratios,
                              deconv=deconv, galalithe=gal)
        assert "Peak Detection" in text
        assert "Diagnostic Ratios" in text
        assert "AMIDE I DECONVOLUTION" in text
        assert "Galalithe" in text

    def test_report_writes_to_file(self, tmp_path, clean_casein_spectrum):
        out = tmp_path / "report.txt"
        text = text_report(clean_casein_spectrum, output=str(out))
        assert out.exists()
        on_disk = out.read_text(encoding="utf-8")
        assert on_disk == text


class TestHtmlReport:

    def test_basic_html(self, clean_casein_spectrum):
        html = html_report(clean_casein_spectrum)
        assert "<!DOCTYPE html>" in html or "<html" in html
        assert "</html>" in html
        assert clean_casein_spectrum.name in html

    def test_html_with_full_sections(self, clean_casein_spectrum,
                                       clean_galalithe_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="all")
        ratios = diagnostic_ratios(clean_casein_spectrum)
        gal = galalithe_indicators(clean_galalithe_spectrum)
        changes = highlight_changes(clean_casein_spectrum,
                                        clean_galalithe_spectrum)
        deconv = deconvolve_amide_I(clean_casein_spectrum)
        html = html_report(clean_casein_spectrum,
                              peaks=peaks, ratios=ratios, deconv=deconv,
                              galalithe=gal, changes=changes)
        assert "Detected peaks" in html
        assert "Diagnostic ratios" in html
        assert "deconvolution" in html.lower()
        assert "Galalithe" in html or "galalithe" in html
        assert "comparison" in html.lower() or "Before / After" in html

    def test_html_writes_to_file(self, tmp_path, clean_casein_spectrum):
        out = tmp_path / "report.html"
        html_report(clean_casein_spectrum, output=str(out))
        assert out.exists()
        content = out.read_text(encoding="utf-8")
        assert "</html>" in content
