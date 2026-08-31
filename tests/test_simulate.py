"""
Tests for casein_ftir.simulate
"""

import numpy as np
import pytest

from casein_ftir.io_module import Spectrum
from casein_ftir.simulate import (
    simulate_casein, simulate_galalithe, gaussian, lorentzian,
    SPECIES_SHIFTS, DEFAULT_FWHM,
)
from casein_ftir.database import CASEIN_BANDS


# ---------------------------------------------------------------------------
# Line shapes
# ---------------------------------------------------------------------------

class TestLineShapes:

    def test_gaussian_peak_at_center(self):
        x = np.linspace(1500, 1800, 1000)
        y = gaussian(x, center=1648.0, amplitude=1.0, fwhm=30.0)
        idx = int(np.argmax(y))
        assert abs(x[idx] - 1648.0) < 1.0

    def test_gaussian_amplitude_at_center(self):
        x = np.array([1648.0])
        y = gaussian(x, 1648.0, amplitude=2.5, fwhm=30.0)
        assert y[0] == pytest.approx(2.5)

    def test_gaussian_fwhm(self):
        """FWHM check : the curve at center ± FWHM/2 must equal half max."""
        x = np.linspace(1500, 1800, 5000)
        y = gaussian(x, 1648.0, 1.0, fwhm=30.0)
        half_max_idx = np.where(y >= 0.5)[0]
        width = x[half_max_idx[-1]] - x[half_max_idx[0]]
        assert width == pytest.approx(30.0, rel=0.02)

    def test_lorentzian_peak_at_center(self):
        x = np.linspace(1500, 1800, 1000)
        y = lorentzian(x, 1648.0, 1.0, fwhm=30.0)
        idx = int(np.argmax(y))
        assert abs(x[idx] - 1648.0) < 1.0


# ---------------------------------------------------------------------------
# Casein simulation
# ---------------------------------------------------------------------------

class TestSimulateCasein:

    def test_returns_spectrum(self, clean_casein_spectrum):
        assert isinstance(clean_casein_spectrum, Spectrum)
        assert clean_casein_spectrum.n_points > 100

    def test_range_respected(self):
        s = simulate_casein(wn_min=800, wn_max=3500, resolution=2.0)
        lo, hi = s.range
        assert lo == pytest.approx(800.0, abs=2)
        assert hi == pytest.approx(3500.0, abs=2)

    def test_amide_I_present(self, clean_casein_spectrum):
        """L'amide I doit dominer la région 1600-1700 cm-1."""
        s = clean_casein_spectrum
        # find max in amide I window
        mask = (s.wavenumber >= 1600) & (s.wavenumber <= 1700)
        peak_wn = s.wavenumber[mask][np.argmax(s.absorbance[mask])]
        assert 1635 <= peak_wn <= 1660

    def test_amide_II_present(self, clean_casein_spectrum):
        s = clean_casein_spectrum
        mask = (s.wavenumber >= 1510) & (s.wavenumber <= 1580)
        peak_wn = s.wavenumber[mask][np.argmax(s.absorbance[mask])]
        assert 1530 <= peak_wn <= 1550

    def test_phosphate_present(self, clean_casein_spectrum):
        """La caséine doit avoir un signal phosphate (~1100 cm-1)."""
        s = clean_casein_spectrum
        mask = (s.wavenumber >= 1050) & (s.wavenumber <= 1150)
        assert s.absorbance[mask].max() > 0.1

    def test_metadata_stored(self, clean_casein_spectrum):
        meta = clean_casein_spectrum.metadata
        assert meta["type"] == "simulated"
        assert meta["species"] == "bovine"
        assert "secondary_structure" in meta

    @pytest.mark.parametrize("species", ["bovine", "caprine", "ovine"])
    def test_species_supported(self, species):
        s = simulate_casein(species=species, noise_level=0.0)
        assert s.metadata["species"] == species

    def test_secondary_structure_sums_to_one_by_default(self):
        s = simulate_casein()
        total = sum(s.metadata["secondary_structure"].values())
        assert total == pytest.approx(1.0, abs=0.001)

    def test_custom_secondary_structure(self):
        struct = {
            "random_coil": 0.10, "alpha_helix": 0.70, "beta_sheet": 0.10,
            "beta_turn": 0.05, "beta_aggregate": 0.025,
            "beta_antiparallel": 0.025,
        }
        s = simulate_casein(secondary_structure=struct, noise_level=0.0)
        assert s.metadata["secondary_structure"] == struct

    def test_noise_increases_variance(self):
        s_clean = simulate_casein(noise_level=0.0)
        s_noisy = simulate_casein(noise_level=0.01)
        # In flat regions (no peak), variance is dominated by noise
        # take 3700-3900 cm-1 zone where there's nothing.
        mask = (s_clean.wavenumber >= 3700) & (s_clean.wavenumber <= 3900)
        assert s_noisy.absorbance[mask].std() > s_clean.absorbance[mask].std()


# ---------------------------------------------------------------------------
# Galalithe simulation
# ---------------------------------------------------------------------------

class TestSimulateGalalithe:

    def test_returns_spectrum(self, clean_galalithe_spectrum):
        assert isinstance(clean_galalithe_spectrum, Spectrum)

    def test_crosslink_increases_CH2_signal(self):
        """Plus on réticule, plus le pic CH2 ~2920 doit être fort."""
        low  = simulate_galalithe(crosslink_degree=0.1, noise_level=0.0)
        high = simulate_galalithe(crosslink_degree=0.9, noise_level=0.0)
        mask = (low.wavenumber >= 2900) & (low.wavenumber <= 2940)
        assert high.absorbance[mask].max() > low.absorbance[mask].max()

    def test_no_crosslink_resembles_casein(self):
        """À taux nul, galalithe ~= caséine (sauf shift de structure 2nd)."""
        s_gal = simulate_galalithe(crosslink_degree=0.0, residual_formol=0.0,
                                     noise_level=0.0)
        s_cas = simulate_casein(noise_level=0.0)
        # spectres alignés
        np.testing.assert_array_equal(s_gal.wavenumber, s_cas.wavenumber)
        # similar overall shape -- correlation > 0.95
        corr = np.corrcoef(s_gal.absorbance, s_cas.absorbance)[0, 1]
        assert corr > 0.95

    def test_residual_formol_adds_1720_peak(self):
        with_f    = simulate_galalithe(crosslink_degree=0.5, residual_formol=0.3,
                                          noise_level=0.0)
        without_f = simulate_galalithe(crosslink_degree=0.5, residual_formol=0.0,
                                          noise_level=0.0)
        mask = (with_f.wavenumber >= 1700) & (with_f.wavenumber <= 1740)
        assert with_f.absorbance[mask].max() > without_f.absorbance[mask].max()
