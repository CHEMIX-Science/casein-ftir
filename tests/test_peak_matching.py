"""
Tests for casein_ftir.reference_module.peak_matching_score
"""

import math
import numpy as np
import pytest

from casein_ftir.io_module import Spectrum
from casein_ftir.simulate import simulate_casein, simulate_galalithe
from casein_ftir.reference_module import (
    peak_matching_score, report_peak_matching, compare_to_reference,
)


# ---------------------------------------------------------------------------
# Basic structure
# ---------------------------------------------------------------------------

class TestPeakMatchingStructure:

    def test_returns_expected_keys(self):
        cas = simulate_casein(noise_level=0.0)
        result = peak_matching_score(cas, cas)
        expected = {"n_ref_peaks", "n_sample_peaks", "n_matched",
                     "match_rate", "weighted_match_rate", "matches",
                     "tolerance_cm-1"}
        assert set(result.keys()) == expected

    def test_match_self_is_complete(self):
        """Matching a spectrum to itself : 100% match."""
        cas = simulate_casein(noise_level=0.0)
        result = peak_matching_score(cas, cas, tolerance=5.0)
        assert result["match_rate"] == pytest.approx(1.0)
        assert result["weighted_match_rate"] == pytest.approx(1.0)
        assert result["n_matched"] == result["n_ref_peaks"]

    def test_each_match_has_required_fields(self):
        cas = simulate_casein(noise_level=0.0)
        result = peak_matching_score(cas, cas)
        for m in result["matches"]:
            assert set(m.keys()) == {
                "ref_wn", "ref_intensity",
                "nearest_sample_wn", "nearest_sample_intensity",
                "delta_wn", "matched",
            }


# ---------------------------------------------------------------------------
# Tolerance behaviour
# ---------------------------------------------------------------------------

class TestPeakMatchingTolerance:

    def test_tight_tolerance_lowers_match_rate(self):
        cas = simulate_casein(noise_level=0.005)
        gal = simulate_galalithe(crosslink_degree=0.5, noise_level=0.005)
        r_tight = peak_matching_score(cas, gal, tolerance=2.0)
        r_loose = peak_matching_score(cas, gal, tolerance=50.0)
        assert r_loose["match_rate"] >= r_tight["match_rate"]

    def test_tolerance_stored_in_output(self):
        cas = simulate_casein(noise_level=0.0)
        r = peak_matching_score(cas, cas, tolerance=12.5)
        assert r["tolerance_cm-1"] == pytest.approx(12.5)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestPeakMatchingEdgeCases:

    def test_flat_reference_yields_zero_peaks(self):
        wn = np.linspace(4000, 600, 200)
        flat = Spectrum(wn, np.zeros_like(wn))
        cas = simulate_casein(noise_level=0.0)
        r = peak_matching_score(cas, flat)
        assert r["n_ref_peaks"] == 0
        assert math.isnan(r["match_rate"])

    def test_synthetic_reference_has_one_of_three_matches(self):
        """An independently constructed peak set has exactly one known match."""
        wn = np.linspace(4000, 600, 1701)
        sample = Spectrum(wn, sum(np.exp(-((wn-c)/10)**2) for c in (1000, 1500)))
        reference = Spectrum(wn, sum(np.exp(-((wn-c)/10)**2) for c in (1000, 2500, 3000)))
        r = peak_matching_score(sample, reference, tolerance=15.0)
        assert r["n_ref_peaks"] == 3
        assert r["n_matched"] == 1
        assert r["match_rate"] == pytest.approx(1/3)


# ---------------------------------------------------------------------------
# Integration with compare_to_reference
# ---------------------------------------------------------------------------

class TestIntegration:

    def test_compare_includes_peak_matching(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.5, noise_level=0.0)
        result = compare_to_reference(cas, gal)
        assert "peak_matching" in result
        pm = result["peak_matching"]
        assert "match_rate" in pm

    def test_report_contains_peak_matching_section(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.5, noise_level=0.0)
        result = compare_to_reference(cas, gal)
        from casein_ftir.reference_module import report_reference_comparison
        text = report_reference_comparison(result)
        assert "PEAK-BY-PEAK MATCHING" in text
        assert "Match rate" in text


# ---------------------------------------------------------------------------
# Pretty-print
# ---------------------------------------------------------------------------

class TestReport:

    def test_report_is_string(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.5, noise_level=0.0)
        r = peak_matching_score(cas, gal)
        text = report_peak_matching(r)
        assert isinstance(text, str)
        assert "PEAK-BY-PEAK" in text
