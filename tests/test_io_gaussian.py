"""
Tests for casein_ftir.io_gaussian
"""

import numpy as np
import pytest
from pathlib import Path

from casein_ftir.io_module import Spectrum
from casein_ftir.io_gaussian import (
    parse_gaussian_log,
    gaussian_log_to_spectrum,
    list_modes,
    get_scale_factor,
    SCALE_FACTORS,
)


# ---------------------------------------------------------------------------
# Fixture : fake Gaussian log
# ---------------------------------------------------------------------------

GAUSSIAN_LOG_CONTENT = """\
 Entering Gaussian System
 ---------------------------------------------------------------
 # B3LYP/6-311+G(2d,p) Opt Freq SCRF=(IEFPCM,Solvent=Water)
 ---------------------------------------------------------------

 Test molecule

 Harmonic frequencies (cm**-1), IR intensities (KM/Mole)
                      1                      2                      3
 Frequencies --     56.4321              123.5612              285.7811
 Red. masses --     1.5234                2.1432                3.5621
 IR Inten    --     0.4521                1.2354                4.5612
                      4                      5                      6
 Frequencies --    412.6543              598.2143              785.4321
 IR Inten    --     8.4321               25.1432               12.5643
                      7                      8                      9
 Frequencies --   1521.4321             1626.7654             1735.8765
 IR Inten    --   285.7654              125.4321              524.6543
 Normal termination
"""


@pytest.fixture
def fake_log(tmp_path):
    p = tmp_path / "fake.log"
    p.write_text(GAUSSIAN_LOG_CONTENT)
    return p


# ---------------------------------------------------------------------------
# Scale factors
# ---------------------------------------------------------------------------

class TestScaleFactors:

    def test_known_levels(self):
        assert get_scale_factor("B3LYP/6-311+G(2d,p)") == pytest.approx(0.9679)
        assert get_scale_factor("B3LYP/6-31G(d)") == pytest.approx(0.9614)

    def test_case_insensitive(self):
        assert get_scale_factor("b3lyp/6-31g(d)") == get_scale_factor("B3LYP/6-31G(d)")

    def test_whitespace_tolerant(self):
        assert get_scale_factor(" B3LYP / 6-31G(d) ") == \
                get_scale_factor("B3LYP/6-31G(d)")

    def test_unknown_returns_one(self):
        assert get_scale_factor("UnknownMethod/whatever") == 1.0


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

class TestParser:

    def test_parses_frequencies(self, fake_log):
        result = parse_gaussian_log(fake_log)
        # 9 modes (3 blocks of 3)
        assert result["n_modes"] == 9
        np.testing.assert_allclose(
            result["frequencies"][:3],
            [56.4321, 123.5612, 285.7811],
        )

    def test_parses_intensities(self, fake_log):
        result = parse_gaussian_log(fake_log)
        np.testing.assert_allclose(
            result["intensities"][:3],
            [0.4521, 1.2354, 4.5612],
        )

    def test_detects_level(self, fake_log):
        result = parse_gaussian_log(fake_log)
        assert "B3LYP" in result["level"]
        assert "6-311" in result["level"]

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            parse_gaussian_log("/no/such/file.log")

    def test_no_frequencies_raises(self, tmp_path):
        p = tmp_path / "empty.log"
        p.write_text("This is not a Gaussian output\nNothing useful here.\n")
        with pytest.raises(ValueError):
            parse_gaussian_log(p)


# ---------------------------------------------------------------------------
# Spectrum generation
# ---------------------------------------------------------------------------

class TestGaussianToSpectrum:

    def test_returns_spectrum(self, fake_log):
        spec = gaussian_log_to_spectrum(fake_log)
        assert isinstance(spec, Spectrum)
        assert spec.n_points > 100

    def test_metadata_filled(self, fake_log):
        spec = gaussian_log_to_spectrum(fake_log)
        meta = spec.metadata
        assert meta["type"] == "dft_simulated"
        assert meta["n_modes_total"] == 9
        assert meta["scale_factor"] == pytest.approx(0.9679, abs=1e-3)

    def test_strongest_mode_peak_position(self, fake_log):
        """Le mode 1735.88 cm-1 (raw) * 0.9679 = 1680.15 cm-1 (scaled)
        doit donner un pic dans le spectre convolué autour de 1680."""
        spec = gaussian_log_to_spectrum(fake_log, fwhm=10.0)
        # 1735.88 has the largest IR intensity (524.6543)
        # find the global maximum in [1650, 1710]
        mask = (spec.wavenumber >= 1650) & (spec.wavenumber <= 1710)
        peak_wn = spec.wavenumber[mask][np.argmax(spec.absorbance[mask])]
        assert 1675.0 <= peak_wn <= 1685.0

    def test_overridden_scale_factor(self, fake_log):
        spec1 = gaussian_log_to_spectrum(fake_log, scale=1.0)
        spec2 = gaussian_log_to_spectrum(fake_log, scale=0.9)
        # With scale=1.0, the strongest peak should be near 1735;
        # with scale=0.9, near 1562.
        mask1 = (spec1.wavenumber >= 1700) & (spec1.wavenumber <= 1770)
        mask2 = (spec2.wavenumber >= 1530) & (spec2.wavenumber <= 1600)
        peak1 = spec1.wavenumber[mask1][np.argmax(spec1.absorbance[mask1])]
        peak2 = spec2.wavenumber[mask2][np.argmax(spec2.absorbance[mask2])]
        assert 1730 <= peak1 <= 1745
        assert 1555 <= peak2 <= 1570

    def test_imaginary_modes_dropped(self, tmp_path):
        """Un mode imaginaire (fréquence négative) doit être ignoré par défaut."""
        content = """\
 # B3LYP/6-31G(d) Opt Freq

 Test
                      1                      2                      3
 Frequencies --   -123.4567             1500.0000             3200.0000
 IR Inten    --    10.0000               50.0000              100.0000
 Normal termination
"""
        p = tmp_path / "imag.log"
        p.write_text(content)
        spec = gaussian_log_to_spectrum(p)
        assert spec.metadata["n_modes_total"] == 3
        assert spec.metadata["n_modes_used"] == 2

    def test_normalised_max_is_one(self, fake_log):
        spec = gaussian_log_to_spectrum(fake_log)
        assert spec.absorbance.max() == pytest.approx(1.0, abs=1e-6)

    def test_lineshape_choice(self, fake_log):
        s_lor = gaussian_log_to_spectrum(fake_log, lineshape="lorentzian",
                                            fwhm=10.0)
        s_gau = gaussian_log_to_spectrum(fake_log, lineshape="gaussian",
                                            fwhm=10.0)
        # different shapes => different spectra (but both peaks at same positions)
        assert not np.allclose(s_lor.absorbance, s_gau.absorbance)


# ---------------------------------------------------------------------------
# list_modes
# ---------------------------------------------------------------------------

class TestListModes:

    def test_returns_string(self, fake_log):
        text = list_modes(fake_log, top_n=5)
        assert isinstance(text, str)
        assert "B3LYP" in text
        assert "Scaled" in text
        # the strongest mode (1735.88, intensity 524.65) should be at rank 1
        # Just check the value 1735 appears (in the harmonic column)
        assert "1735" in text
