"""
Tests for casein_ftir.quantification
"""

import numpy as np
import pytest

from casein_ftir.io_module import Spectrum
from casein_ftir.quantification import (
    integrate_band, integrate_all_bands, diagnostic_ratios,
    galalithe_indicators,
    beer_lambert_concentration, build_calibration_curve, apply_calibration,
)
from casein_ftir.simulate import simulate_casein


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------

class TestIntegration:

    def test_integral_of_gaussian_is_analytical(self):
        """Integrale d'une gaussienne d'aire connue = sigma*sqrt(2*pi)."""
        wn = np.linspace(1700, 1600, 5000)  # descending
        center = 1648.0
        amp = 1.0
        fwhm = 30.0
        sigma = fwhm / 2.3548
        ab = amp * np.exp(-((wn - center) / sigma) ** 2 / 2.0)
        spec = Spectrum(wn, ab)
        area = integrate_band(spec, 1600, 1700, baseline="none")
        expected = amp * sigma * np.sqrt(2 * np.pi)
        assert area == pytest.approx(expected, rel=0.01)

    def test_linear_baseline_subtraction(self):
        """Test : si on ajoute une ligne de base constante, on retrouve
        l'aire du pic seul après correction linéaire."""
        wn = np.linspace(1700, 1600, 2000)
        amp = 1.0; fwhm = 30.0
        sigma = fwhm / 2.3548
        peak = amp * np.exp(-((wn - 1648) / sigma) ** 2 / 2.0)
        background = 0.3
        spec_with_bg = Spectrum(wn, peak + background)
        spec_only    = Spectrum(wn, peak)
        a_corr = integrate_band(spec_with_bg, 1600, 1700, baseline="linear")
        a_only = integrate_band(spec_only,    1600, 1700, baseline="none")
        assert a_corr == pytest.approx(a_only, rel=0.02)

    def test_empty_window_returns_zero(self, clean_casein_spectrum):
        # Window entirely outside spectrum range
        a = integrate_band(clean_casein_spectrum, 10000, 11000)
        assert a == 0.0

    def test_integrate_all_bands_returns_dict(self, clean_casein_spectrum):
        d = integrate_all_bands(clean_casein_spectrum, band_set="casein")
        assert "amide_I" in d
        assert "amide_II" in d
        assert d["amide_I"] > 0
        assert d["amide_II"] > 0


# ---------------------------------------------------------------------------
# Diagnostic ratios
# ---------------------------------------------------------------------------

class TestRatios:

    def test_ratios_keys_present(self, clean_casein_spectrum):
        r = diagnostic_ratios(clean_casein_spectrum)
        expected_keys = {"A_II_over_I", "A_phos_over_I", "A_CH2_over_I",
                          "A_amideA_over_I", "A_COO_over_I"}
        assert set(r.keys()) == expected_keys

    def test_ratios_are_positive(self, clean_casein_spectrum):
        r = diagnostic_ratios(clean_casein_spectrum)
        for k, v in r.items():
            assert v > 0, f"Ratio {k} should be positive, got {v}"

    def test_A_II_over_I_in_expected_range(self, clean_casein_spectrum):
        """Pour la caséine native, A_II/A_I doit être 0.4 - 1.0 environ."""
        r = diagnostic_ratios(clean_casein_spectrum)
        assert 0.3 < r["A_II_over_I"] < 1.2

    def test_galalithe_indicators_increase_with_crosslink(self):
        from casein_ftir.simulate import simulate_galalithe
        s_low  = simulate_galalithe(crosslink_degree=0.1, noise_level=0.0)
        s_high = simulate_galalithe(crosslink_degree=0.9, noise_level=0.0)
        ind_low  = galalithe_indicators(s_low)
        ind_high = galalithe_indicators(s_high)
        # Crosslink index must strictly grow with crosslink degree
        assert ind_high["crosslink_index"] > ind_low["crosslink_index"]


# ---------------------------------------------------------------------------
# Beer-Lambert
# ---------------------------------------------------------------------------

class TestBeerLambert:

    def test_beer_lambert_basic(self, clean_casein_spectrum):
        # epsilon = 1 (arbitrary), l = 1 cm -> c = A
        c = beer_lambert_concentration(
            clean_casein_spectrum, band_name="amide_II",
            molar_absorptivity=1.0, path_length=1.0, method="peak",
        )
        assert c > 0.0

    def test_beer_lambert_inverse_scaling(self, clean_casein_spectrum):
        """c doit doubler quand epsilon est divisé par 2."""
        c1 = beer_lambert_concentration(
            clean_casein_spectrum, "amide_II", molar_absorptivity=2.0,
            path_length=1.0, method="area",
        )
        c2 = beer_lambert_concentration(
            clean_casein_spectrum, "amide_II", molar_absorptivity=1.0,
            path_length=1.0, method="area",
        )
        assert c2 == pytest.approx(2 * c1, rel=1e-6)

    def test_unknown_band_raises(self, clean_casein_spectrum):
        with pytest.raises(KeyError):
            beer_lambert_concentration(
                clean_casein_spectrum, "not_a_band", 1.0, 1.0
            )

    def test_zero_path_returns_nan(self, clean_casein_spectrum):
        c = beer_lambert_concentration(
            clean_casein_spectrum, "amide_II",
            molar_absorptivity=1.0, path_length=0.0
        )
        assert np.isnan(c)


# ---------------------------------------------------------------------------
# Calibration curve
# ---------------------------------------------------------------------------

class TestCalibration:

    def test_linear_calibration_recovers_slope(self):
        """On simule une série de spectres caséine avec amplitude variable
        et on vérifie que la régression linéaire retrouve cette amplitude."""
        # Build a series of caseins with different overall intensities
        # by adding a constant scalar in absorbance.
        base = simulate_casein(noise_level=0.0)
        concs = [1.0, 2.0, 3.0, 4.0, 5.0]
        spectra = []
        for c in concs:
            s = base.copy()
            s.absorbance = s.absorbance * c
            spectra.append(s)
        cal = build_calibration_curve(spectra, concs, band_name="amide_II",
                                          method="area")
        # R^2 must be near 1, intercept near zero, slope > 0
        assert cal["r_squared"] > 0.999
        assert cal["slope"] > 0
        assert abs(cal["intercept"]) < 1.0

    def test_apply_calibration_roundtrip(self):
        """Après calibration, on doit retrouver la concentration injectée."""
        base = simulate_casein(noise_level=0.0)
        concs = [1.0, 2.0, 3.0, 4.0, 5.0]
        spectra = []
        for c in concs:
            s = base.copy()
            s.absorbance = s.absorbance * c
            spectra.append(s)
        cal = build_calibration_curve(spectra, concs, band_name="amide_II",
                                          method="area")
        # Verify : the predicted concentration must equal the input
        for spec, c_true in zip(spectra, concs):
            c_pred = apply_calibration(spec, cal)
            assert c_pred == pytest.approx(c_true, rel=0.02)

    def test_calibration_mismatched_lengths_raises(self):
        with pytest.raises(ValueError):
            build_calibration_curve([], [1.0, 2.0])
