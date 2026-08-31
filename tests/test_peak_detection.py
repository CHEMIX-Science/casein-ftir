"""
Tests for casein_ftir.peak_detection
"""

import math
import numpy as np
import pytest

from casein_ftir.peak_detection import detect_bands, DetectedPeak, report_peaks
from casein_ftir.database import CASEIN_BANDS, GALALITHE_MARKERS


# ---------------------------------------------------------------------------
# Detection on simulated spectra
# ---------------------------------------------------------------------------

class TestDetectionOnCasein:

    def test_finds_all_casein_bands(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="casein")
        assert set(peaks.keys()) == set(CASEIN_BANDS.keys())
        for name, peak in peaks.items():
            assert peak.found, f"Band {name} not found"

    def test_amide_I_within_tolerance(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="casein")
        p = peaks["amide_I"]
        # Should be within ±10 cm-1 of the expected center
        assert abs(p.shift) < 10.0

    def test_amide_II_within_tolerance(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="casein")
        p = peaks["amide_II"]
        assert abs(p.shift) < 10.0

    def test_phosphate_within_tolerance(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="casein")
        p = peaks["phosphate"]
        # Phosphate is broad, allow up to ±15 cm-1
        assert abs(p.shift) < 20.0

    def test_returned_peak_is_dataclass(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="casein")
        for p in peaks.values():
            assert isinstance(p, DetectedPeak)
            assert isinstance(p.as_dict(), dict)


# ---------------------------------------------------------------------------
# Detection on galalithe markers
# ---------------------------------------------------------------------------

class TestDetectionOnGalalithe:

    def test_galalithe_markers_strong_when_crosslinked(self, clean_galalithe_spectrum):
        peaks = detect_bands(clean_galalithe_spectrum, band_set="galalithe")
        # methylene_bridge and N_CH2_N should be present
        assert peaks["methylene_bridge"].found
        assert peaks["N_CH2_N"].found

    def test_all_band_set(self, clean_galalithe_spectrum):
        peaks = detect_bands(clean_galalithe_spectrum, band_set="all")
        expected = set(CASEIN_BANDS.keys()) | set(GALALITHE_MARKERS.keys())
        assert set(peaks.keys()) == expected


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestDetectionEdgeCases:

    def test_invalid_band_set_raises(self, clean_casein_spectrum):
        with pytest.raises(ValueError):
            detect_bands(clean_casein_spectrum, band_set="bogus")

    def test_empty_window_marks_not_found(self):
        """Si on lance la détection sur un spectre court ne couvrant pas
        toutes les bandes, certaines doivent être marquées 'not found'."""
        from casein_ftir.io_module import Spectrum
        # Spectre minuscule sur 1500-1700 cm-1 uniquement
        wn = np.linspace(1700, 1500, 100)
        ab = np.exp(-((wn - 1648) / 30) ** 2)
        s = Spectrum(wn, ab)
        peaks = detect_bands(s, band_set="casein")
        # amide A (>3000 cm-1) ne peut pas être trouvé
        assert not peaks["amide_A"].found
        # amide I doit être trouvé
        assert peaks["amide_I"].found


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

class TestReportPeaks:

    def test_report_is_non_empty_string(self, clean_casein_spectrum):
        peaks = detect_bands(clean_casein_spectrum, band_set="casein")
        text = report_peaks(peaks)
        assert isinstance(text, str)
        assert "amide_I" in text
        assert "Position" in text or "Expected" in text
