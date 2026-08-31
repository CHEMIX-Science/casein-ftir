"""
Tests for casein_ftir.reference_module
"""

import numpy as np
import pytest

from casein_ftir.io_module import save_spectrum_csv, load_spectrum
from casein_ftir.simulate import simulate_casein, simulate_galalithe
from casein_ftir.reference_module import (
    load_reference,
    compare_to_reference,
    subtract_reference,
    report_reference_comparison,
)


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

class TestLoadReference:

    def test_load_with_preprocess(self, tmp_path, clean_casein_spectrum):
        path = tmp_path / "ref.csv"
        save_spectrum_csv(clean_casein_spectrum, path)
        ref = load_reference(str(path), preprocess=True)
        assert ref.metadata.get("is_reference") is True
        assert "normalization" in ref.metadata
        assert "baseline_method" in ref.metadata

    def test_load_without_preprocess(self, tmp_path, clean_casein_spectrum):
        path = tmp_path / "ref.csv"
        save_spectrum_csv(clean_casein_spectrum, path)
        ref = load_reference(str(path), preprocess=False)
        assert ref.metadata.get("is_reference") is True
        # should NOT have preprocess metadata flags
        assert "normalization" not in ref.metadata


# ---------------------------------------------------------------------------
# Compare to reference
# ---------------------------------------------------------------------------

class TestCompareToReference:

    def test_identical_spectra_high_quality(self, clean_casein_spectrum):
        result = compare_to_reference(clean_casein_spectrum,
                                          clean_casein_spectrum)
        assert result["overall_quality"] == "excellent"
        assert result["similarity"]["pearson_corr"] == pytest.approx(1.0)

    def test_galalithe_vs_casein_band_level_changes(self):
        """Globalement, caséine et galalithe restent très corrélées
        (>0.95 Pearson) parce que le squelette protéique est conservé.
        Mais des bandes spécifiques doivent montrer des excès."""
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.9, residual_formol=0.3,
                                    noise_level=0.0)
        result = compare_to_reference(gal, cas)
        # The CH2 / methylene bridge should be in excess
        flags = {n: d["flag"] for n, d in result["band_deviations"].items()}
        # at least one of the galalithe-specific bands flagged as excess
        gal_flags = [flags["methylene_bridge"], flags["free_aldehyde_CO"],
                       flags["N_CH2_N"]]
        assert any(f in ("excess", "sample-only") for f in gal_flags)

    def test_per_band_flags_make_sense(self):
        cas = simulate_casein(noise_level=0.0)
        gal = simulate_galalithe(crosslink_degree=0.9, residual_formol=0.0,
                                    noise_level=0.0)
        result = compare_to_reference(gal, cas, band_change_threshold=0.10)
        # methylene_bridge should be flagged as 'excess'
        flag = result["band_deviations"]["methylene_bridge"]["flag"]
        assert flag in ("excess", "sample-only")

    def test_result_structure(self, clean_casein_spectrum):
        result = compare_to_reference(clean_casein_spectrum,
                                          clean_casein_spectrum)
        assert set(result.keys()) == {
            "similarity", "band_deviations",
            "overall_quality", "recommended_action",
            "peak_matching",
        }
        for name, d in result["band_deviations"].items():
            assert "sample_area" in d
            assert "reference_area" in d
            assert "rel_area_diff_pct" in d
            assert "shift_cm-1" in d
            assert "flag" in d


# ---------------------------------------------------------------------------
# Subtract reference
# ---------------------------------------------------------------------------

class TestSubtractReference:

    def test_subtract_self_is_zero(self, clean_casein_spectrum):
        diff = subtract_reference(clean_casein_spectrum,
                                       clean_casein_spectrum, scale=1.0)
        assert np.allclose(diff.absorbance, 0.0, atol=1e-9)

    def test_subtract_with_scale(self, clean_casein_spectrum):
        diff = subtract_reference(clean_casein_spectrum,
                                       clean_casein_spectrum, scale=0.5)
        # we should get half the spectrum back
        np.testing.assert_allclose(
            diff.absorbance, 0.5 * clean_casein_spectrum.absorbance,
            rtol=1e-6,
        )

    def test_subtract_metadata(self, clean_casein_spectrum):
        diff = subtract_reference(clean_casein_spectrum,
                                       clean_casein_spectrum, scale=1.0)
        assert diff.metadata["type"] == "subtracted"
        assert diff.metadata["scale_factor"] == 1.0


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

class TestReportReferenceComparison:

    def test_report_contains_main_sections(self, clean_casein_spectrum):
        result = compare_to_reference(clean_casein_spectrum,
                                          clean_casein_spectrum)
        text = report_reference_comparison(result)
        assert "REFERENCE COMPARISON" in text
        assert "Pearson" in text
        assert "Per-band deviations" in text
        assert "amide_I" in text

    def test_report_is_string(self, clean_casein_spectrum):
        result = compare_to_reference(clean_casein_spectrum,
                                          clean_casein_spectrum)
        assert isinstance(report_reference_comparison(result), str)
