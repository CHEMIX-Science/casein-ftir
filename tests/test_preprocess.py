"""
Tests for casein_ftir.preprocess
"""

import numpy as np
import pytest

from casein_ftir.io_module import Spectrum
from casein_ftir.preprocess import (
    baseline_als, baseline_polynomial, baseline_rubberband,
    smooth_savgol,
    normalize_minmax, normalize_vector, normalize_area, normalize_to_band,
    full_preprocess,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _peak_position(spec, wn_min, wn_max):
    mask = (spec.wavenumber >= wn_min) & (spec.wavenumber <= wn_max)
    return spec.wavenumber[mask][np.argmax(spec.absorbance[mask])]


# ---------------------------------------------------------------------------
# Baseline correction
# ---------------------------------------------------------------------------

class TestBaseline:

    def test_polynomial_removes_linear_drift(self, clean_casein_spectrum):
        """Une ligne de base polynomiale doit aplatir une dérive linéaire."""
        s = clean_casein_spectrum.copy()
        drift = 0.5 * (s.wavenumber - s.wavenumber.min()) / np.ptp(s.wavenumber)
        s.absorbance = s.absorbance + drift
        s_corr = baseline_polynomial(s, degree=2)
        # endpoints should now be similar
        end_diff = abs(s_corr.absorbance[0] - s_corr.absorbance[-1])
        assert end_diff < 0.2

    def test_als_does_not_destroy_peaks(self, clean_casein_spectrum):
        """ALS doit garder les pics intacts."""
        peak_before = _peak_position(clean_casein_spectrum, 1600, 1700)
        s_corr = baseline_als(clean_casein_spectrum, lam=1e5, p=0.001)
        peak_after = _peak_position(s_corr, 1600, 1700)
        assert abs(peak_before - peak_after) < 5.0

    def test_rubberband_returns_non_negative_for_positive_input(self):
        wn = np.linspace(4000, 600, 500)
        ab = np.exp(-((wn - 1648) / 50) ** 2) + 0.05
        s = Spectrum(wn, ab)
        s_corr = baseline_rubberband(s)
        # After rubberband, baseline endpoints should be near zero
        assert s_corr.absorbance.min() >= -0.1

    def test_baseline_records_metadata(self, clean_casein_spectrum):
        s = baseline_polynomial(clean_casein_spectrum, degree=3)
        assert "baseline_method" in s.metadata


# ---------------------------------------------------------------------------
# Smoothing
# ---------------------------------------------------------------------------

class TestSmoothing:

    def test_smooth_reduces_noise_variance(self, noisy_casein_spectrum):
        s_smooth = smooth_savgol(noisy_casein_spectrum, window=15, poly=3)
        # check noise reduction in a flat region (3700-3900 cm-1)
        mask = (noisy_casein_spectrum.wavenumber >= 3700) & \
               (noisy_casein_spectrum.wavenumber <= 3900)
        std_before = noisy_casein_spectrum.absorbance[mask].std()
        std_after  = s_smooth.absorbance[mask].std()
        assert std_after < std_before

    def test_smooth_preserves_peak_position(self, noisy_casein_spectrum):
        peak_before = _peak_position(noisy_casein_spectrum, 1600, 1700)
        s_smooth = smooth_savgol(noisy_casein_spectrum, window=11, poly=3)
        peak_after = _peak_position(s_smooth, 1600, 1700)
        assert abs(peak_before - peak_after) < 4.0

    def test_smooth_handles_even_window(self, clean_casein_spectrum):
        """Si on passe un window pair, ça doit être corrigé automatiquement."""
        s = smooth_savgol(clean_casein_spectrum, window=10, poly=3)
        # ne doit pas planter
        assert s.n_points == clean_casein_spectrum.n_points


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

class TestNormalization:

    def test_minmax(self, clean_casein_spectrum):
        s = normalize_minmax(clean_casein_spectrum)
        assert s.absorbance.min() == pytest.approx(0.0, abs=1e-6)
        assert s.absorbance.max() == pytest.approx(1.0, abs=1e-6)

    def test_vector_norm(self, clean_casein_spectrum):
        s = normalize_vector(clean_casein_spectrum)
        assert np.linalg.norm(s.absorbance) == pytest.approx(1.0, abs=1e-6)

    def test_area_normalisation(self, clean_casein_spectrum):
        s = normalize_area(clean_casein_spectrum)
        area = float(np.trapezoid(s.absorbance, -s.wavenumber))
        assert area == pytest.approx(1.0, abs=1e-6)

    def test_normalize_to_band_sets_peak_to_one(self, clean_casein_spectrum):
        s = normalize_to_band(clean_casein_spectrum, center=1648, halfwidth=30)
        mask = (s.wavenumber >= 1618) & (s.wavenumber <= 1678)
        assert s.absorbance[mask].max() == pytest.approx(1.0, abs=1e-6)

    def test_normalize_to_band_outside_range_raises(self):
        wn = np.linspace(2000, 1000, 100)
        ab = np.ones(100)
        s = Spectrum(wn, ab)
        with pytest.raises(ValueError):
            normalize_to_band(s, center=3500, halfwidth=20)


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------

class TestFullPreprocess:

    def test_default_pipeline_runs(self, noisy_casein_spectrum):
        s = full_preprocess(noisy_casein_spectrum)
        assert s.n_points == noisy_casein_spectrum.n_points
        # amide I has been normalised to 1
        peak_pos = _peak_position(s, 1600, 1700)
        peak_val = s.absorbance[(s.wavenumber >= 1618) & (s.wavenumber <= 1678)].max()
        assert peak_val == pytest.approx(1.0, abs=0.05)

    def test_pipeline_with_no_baseline(self, noisy_casein_spectrum):
        s = full_preprocess(noisy_casein_spectrum, baseline=None,
                              smooth=True, normalize=None)
        assert s.n_points == noisy_casein_spectrum.n_points

    def test_invalid_baseline_raises(self, clean_casein_spectrum):
        with pytest.raises(ValueError):
            full_preprocess(clean_casein_spectrum, baseline="not_a_method")

    def test_invalid_normalize_raises(self, clean_casein_spectrum):
        with pytest.raises(ValueError):
            full_preprocess(clean_casein_spectrum, normalize="nope")
