"""
Tests for casein_ftir.amide_deconv

Stratégie centrale : on simule un spectre avec une structure secondaire
connue, on déconvolue, on vérifie qu'on retrouve à peu près les fractions
injectées.
"""

import numpy as np
import pytest

from casein_ftir.simulate import simulate_casein
from casein_ftir.amide_deconv import (deconvolve_amide_I, report_deconvolution,
                                         _gaussian)
from casein_ftir.database import AMIDE_I_SUBBANDS


# ---------------------------------------------------------------------------
# Output shape and keys
# ---------------------------------------------------------------------------

class TestDeconvolutionOutput:

    def test_returns_expected_keys(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        expected_keys = {
            "components", "structure_pct_per_subband", "summary",
            "total_area", "fit_quality_r2", "residual",
            "x", "y_corrected", "y_fit",
        }
        assert set(result.keys()) == expected_keys

    def test_all_subbands_present(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        assert set(result["components"].keys()) == set(AMIDE_I_SUBBANDS.keys())

    def test_each_component_has_required_fields(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        for name, comp in result["components"].items():
            assert "center" in comp
            assert "fwhm" in comp
            assert "area" in comp
            assert comp["area"] >= 0
            assert comp["fwhm"] > 0

    def test_summary_lumped_categories(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        s = result["summary"]
        for key in ("alpha_helix_pct", "beta_sheet_pct", "random_coil_pct",
                     "beta_turn_pct"):
            assert key in s

    def test_percentages_sum_to_100(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        total = sum(result["structure_pct_per_subband"].values())
        assert total == pytest.approx(100.0, abs=0.5)


# ---------------------------------------------------------------------------
# Fit quality
# ---------------------------------------------------------------------------

class TestFitQuality:

    def test_R2_near_1_on_clean_spectrum(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        assert result["fit_quality_r2"] > 0.99

    def test_R2_still_good_on_noisy_spectrum(self, noisy_casein_spectrum):
        result = deconvolve_amide_I(noisy_casein_spectrum)
        # noisy 0.005 stddev, still high quality fit
        assert result["fit_quality_r2"] > 0.95


# ---------------------------------------------------------------------------
# Recovery of injected secondary structure
# ---------------------------------------------------------------------------

class TestSecondaryStructureRecovery:
    """On vérifie qu'on retrouve la structure injectée dans la simulation."""

    def test_random_coil_dominant_in_native_casein(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        pcts = result["structure_pct_per_subband"]
        # Random coil should be the largest contribution
        max_name = max(pcts, key=pcts.get)
        assert max_name == "random_coil"

    def test_dominantly_alpha_recovered(self):
        """Si on simule une protéine 70% alpha, on doit retrouver
        ~70% alpha à la déconvolution."""
        struct = {
            "random_coil": 0.10, "alpha_helix": 0.70, "beta_sheet": 0.05,
            "beta_turn": 0.05, "beta_aggregate": 0.05,
            "beta_antiparallel": 0.05,
        }
        s = simulate_casein(secondary_structure=struct, noise_level=0.0)
        result = deconvolve_amide_I(s)
        alpha_pct = result["structure_pct_per_subband"]["alpha_helix"]
        # tolérance large : 70% injecté -> entre 55% et 85%
        assert 55.0 < alpha_pct < 85.0

    def test_dominantly_beta_recovered(self):
        struct = {
            "random_coil": 0.10, "alpha_helix": 0.10, "beta_sheet": 0.55,
            "beta_turn": 0.05, "beta_aggregate": 0.10,
            "beta_antiparallel": 0.10,
        }
        s = simulate_casein(secondary_structure=struct, noise_level=0.0)
        result = deconvolve_amide_I(s)
        beta = result["summary"]["beta_sheet_pct"]
        assert beta > 50.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class TestReport:

    def test_report_deconvolution_returns_string(self, clean_casein_spectrum):
        result = deconvolve_amide_I(clean_casein_spectrum)
        text = report_deconvolution(result)
        assert isinstance(text, str)
        assert "alpha_helix" in text
        assert "R^2" in text or "R2" in text or "R²" in text or "Fit" in text


class TestInternalGaussian:

    def test_gaussian_at_center(self):
        x = np.array([1648.0])
        y = _gaussian(x, c=1648.0, a=2.0, fwhm=30.0)
        assert y[0] == pytest.approx(2.0)


class TestExportCsv:

    def test_export_creates_file(self, tmp_path, clean_casein_spectrum):
        from casein_ftir.amide_deconv import (deconvolve_amide_I,
                                                  export_deconvolution_csv)
        result = deconvolve_amide_I(clean_casein_spectrum)
        out = tmp_path / "deconv.csv"
        export_deconvolution_csv(result, str(out))
        assert out.exists()

    def test_exported_csv_has_correct_columns(self, tmp_path, clean_casein_spectrum):
        from casein_ftir.amide_deconv import (deconvolve_amide_I,
                                                  export_deconvolution_csv)
        result = deconvolve_amide_I(clean_casein_spectrum)
        out = tmp_path / "deconv.csv"
        export_deconvolution_csv(result, str(out))
        # parse the header
        content = out.read_text()
        header_line = [l for l in content.splitlines()
                        if not l.startswith("#")][0]
        cols = header_line.split(",")
        assert "wavenumber" in cols
        assert "observed" in cols
        assert "fit_total" in cols
        assert "residual" in cols
        # all 6 sub-bands present
        for sb in ("beta_aggregate", "beta_sheet", "random_coil",
                    "alpha_helix", "beta_turn", "beta_antiparallel"):
            assert sb in cols

    def test_sum_of_subbands_equals_fit(self, tmp_path, clean_casein_spectrum):
        """For each x, the sum of the 6 sub-band columns should equal
        the fit_total column."""
        from casein_ftir.amide_deconv import (deconvolve_amide_I,
                                                  export_deconvolution_csv)
        result = deconvolve_amide_I(clean_casein_spectrum)
        out = tmp_path / "deconv.csv"
        export_deconvolution_csv(result, str(out))
        # parse
        rows = []
        with out.open() as f:
            for line in f:
                if line.startswith("#"):
                    continue
                rows.append(line.strip().split(","))
        header = rows[0]
        data = np.array([[float(v) for v in r] for r in rows[1:]])
        sub_band_cols = [header.index(sb) for sb in
                          ("beta_aggregate", "beta_sheet", "random_coil",
                           "alpha_helix", "beta_turn", "beta_antiparallel")]
        fit_col = header.index("fit_total")
        sum_subs = data[:, sub_band_cols].sum(axis=1)
        fit = data[:, fit_col]
        np.testing.assert_allclose(sum_subs, fit, rtol=1e-3)
