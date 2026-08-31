"""
Tests for casein_ftir.compare
"""

import numpy as np
import pytest

from casein_ftir.io_module import Spectrum
from casein_ftir.simulate import simulate_casein, simulate_galalithe
from casein_ftir.compare import (
    resample_to_common_axis, difference_spectrum, highlight_changes,
    similarity_metric, report_changes,
)


# ---------------------------------------------------------------------------
# Resampling
# ---------------------------------------------------------------------------

class TestResample:

    def test_resample_aligns_axes(self):
        s1 = simulate_casein(wn_min=600, wn_max=4000, resolution=2.0,
                                noise_level=0.0)
        s2 = simulate_casein(wn_min=800, wn_max=3800, resolution=4.0,
                                noise_level=0.0)
        a, b = resample_to_common_axis(s1, s2, resolution=2.0)
        np.testing.assert_array_equal(a.wavenumber, b.wavenumber)

    def test_resample_descending(self):
        s1 = simulate_casein(noise_level=0.0)
        s2 = simulate_galalithe(crosslink_degree=0.5, noise_level=0.0)
        a, _ = resample_to_common_axis(s1, s2)
        assert a.wavenumber[0] > a.wavenumber[-1]


# ---------------------------------------------------------------------------
# Difference spectrum
# ---------------------------------------------------------------------------

class TestDifferenceSpectrum:

    def test_identical_spectra_give_zero_difference(self, clean_casein_spectrum):
        diff = difference_spectrum(clean_casein_spectrum,
                                       clean_casein_spectrum, normalize=False)
        assert np.allclose(diff.absorbance, 0.0, atol=1e-9)

    def test_galalithe_minus_casein_positive_at_CH2(self):
        """Le pic CH2 ~2920 doit être positif dans la différence
        galalithe - caséine."""
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.7, residual_formol=0.0,
                                    noise_level=0.0)
        diff = difference_spectrum(cas, gal, normalize=False)
        mask = (diff.wavenumber >= 2900) & (diff.wavenumber <= 2940)
        assert diff.absorbance[mask].max() > 0


# ---------------------------------------------------------------------------
# Highlight changes
# ---------------------------------------------------------------------------

class TestHighlightChanges:

    def test_no_change_when_same_spectrum(self, clean_casein_spectrum):
        changes = highlight_changes(clean_casein_spectrum,
                                        clean_casein_spectrum, threshold=0.05)
        for name, c in changes.items():
            # rel_change should be very small in all bands present
            if not np.isnan(c["rel_change"]) and c["rel_change"] != float("inf"):
                assert abs(c["rel_change"]) < 0.05

    def test_galalithe_markers_increase_after_crosslink(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.8, residual_formol=0.0,
                                    noise_level=0.0)
        changes = highlight_changes(cas, gal, threshold=0.10)
        # methylene_bridge band should be flagged as 'increased' or 'APPEARED'
        flag = changes["methylene_bridge"]["flag"]
        assert flag in ("increased", "APPEARED")

    def test_residual_formol_appears(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.5, residual_formol=0.3,
                                    noise_level=0.0)
        changes = highlight_changes(cas, gal)
        flag = changes["free_aldehyde_CO"]["flag"]
        assert flag in ("increased", "APPEARED")


# ---------------------------------------------------------------------------
# Similarity metrics
# ---------------------------------------------------------------------------

class TestSimilarity:

    def test_self_similarity_is_one(self, clean_casein_spectrum):
        sim = similarity_metric(clean_casein_spectrum, clean_casein_spectrum)
        assert sim["pearson_corr"] == pytest.approx(1.0, abs=1e-6)
        assert sim["cosine"] == pytest.approx(1.0, abs=1e-6)
        assert sim["rmse"] == pytest.approx(0.0, abs=1e-9)
        assert sim["spectral_overlap"] == pytest.approx(1.0, abs=1e-6)

    def test_casein_galalithe_still_similar(self):
        """Caséine et galalithe partagent la même backbone : Pearson > 0.9."""
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.5, residual_formol=0.05,
                                    noise_level=0.0)
        sim = similarity_metric(cas, gal)
        assert sim["pearson_corr"] > 0.9

    def test_metrics_keys(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.5, noise_level=0.0)
        sim = similarity_metric(cas, gal)
        expected = {"pearson_corr", "cosine", "rmse", "spectral_overlap"}
        assert set(sim.keys()) == expected


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

class TestReportChanges:

    def test_report_changes_returns_string(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.5, noise_level=0.0)
        changes = highlight_changes(cas, gal)
        text = report_changes(changes)
        assert isinstance(text, str)
        assert "methylene_bridge" in text
        assert len(text.splitlines()) > 5
