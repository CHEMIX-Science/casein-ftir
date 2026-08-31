"""
Tests for casein_ftir.io_module
"""

import numpy as np
import pytest

from casein_ftir.io_module import Spectrum, load_spectrum, save_spectrum_csv


# ---------------------------------------------------------------------------
# Spectrum class
# ---------------------------------------------------------------------------

class TestSpectrum:

    def test_construct_with_arrays(self):
        wn = np.array([3000.0, 2000.0, 1000.0])
        ab = np.array([0.1, 0.5, 0.3])
        s = Spectrum(wn, ab, name="test")
        assert s.n_points == 3
        assert s.name == "test"

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError):
            Spectrum(np.array([1, 2, 3]), np.array([0.1, 0.2]))

    def test_ascending_input_gets_reordered_descending(self):
        """Convention IR : axe descendant."""
        wn = np.array([1000.0, 2000.0, 3000.0])   # ascending
        ab = np.array([0.3, 0.5, 0.1])
        s = Spectrum(wn, ab)
        # After init, wavenumber must be descending
        assert s.wavenumber[0] > s.wavenumber[-1]
        # And absorbance must have been flipped accordingly
        assert s.absorbance[0] == pytest.approx(0.1)
        assert s.absorbance[-1] == pytest.approx(0.3)

    def test_descending_input_preserved(self):
        wn = np.array([3000.0, 2000.0, 1000.0])
        ab = np.array([0.1, 0.5, 0.3])
        s = Spectrum(wn, ab)
        np.testing.assert_array_equal(s.wavenumber, wn)
        np.testing.assert_array_equal(s.absorbance, ab)

    def test_range_property(self):
        wn = np.linspace(4000, 600, 100)
        ab = np.zeros_like(wn)
        s = Spectrum(wn, ab)
        lo, hi = s.range
        assert lo == pytest.approx(600.0)
        assert hi == pytest.approx(4000.0)

    def test_slice(self):
        wn = np.linspace(4000, 600, 1000)
        ab = wn / 4000.0
        s = Spectrum(wn, ab)
        sub = s.slice(1500.0, 1800.0)
        assert sub.wavenumber.min() >= 1500.0
        assert sub.wavenumber.max() <= 1800.0
        assert sub.n_points < s.n_points

    def test_copy_is_deep(self):
        wn = np.array([3000.0, 2000.0, 1000.0])
        ab = np.array([0.1, 0.5, 0.3])
        s = Spectrum(wn, ab)
        s2 = s.copy()
        s2.absorbance[0] = 99.0
        assert s.absorbance[0] == pytest.approx(0.1)

    def test_repr_does_not_crash(self):
        wn = np.array([3000.0, 2000.0, 1000.0])
        ab = np.array([0.1, 0.5, 0.3])
        s = Spectrum(wn, ab, name="foo")
        r = repr(s)
        assert "foo" in r
        assert "3000" in r


# ---------------------------------------------------------------------------
# CSV reading
# ---------------------------------------------------------------------------

class TestCsvIO:

    def test_load_simple_csv(self, temp_csv):
        s = load_spectrum(str(temp_csv))
        assert isinstance(s, Spectrum)
        assert s.n_points > 100
        assert s.range[1] > 1000  # max wavenumber reasonably large

    def test_roundtrip_csv(self, tmp_path, clean_casein_spectrum):
        out = tmp_path / "rt.csv"
        save_spectrum_csv(clean_casein_spectrum, out)
        loaded = load_spectrum(str(out))
        assert loaded.n_points == clean_casein_spectrum.n_points
        np.testing.assert_allclose(
            loaded.wavenumber, clean_casein_spectrum.wavenumber, rtol=1e-3
        )
        np.testing.assert_allclose(
            loaded.absorbance, clean_casein_spectrum.absorbance, atol=1e-5
        )

    def test_load_unknown_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_spectrum(str(tmp_path / "nope.csv"))

    def test_tab_separated_csv(self, tmp_path):
        path = tmp_path / "tabbed.tsv"
        with open(path, "w") as f:
            f.write("wn\tab\n")
            for w, a in [(3000, 0.1), (2000, 0.5), (1500, 0.3)]:
                f.write(f"{w}\t{a}\n")
        s = load_spectrum(str(path))
        assert s.n_points == 3
