"""
Shared pytest fixtures for the casein_ftir test suite.
"""

import os
import sys

# Make sure the package is importable even if not pip-installed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np
import pytest

from casein_ftir.io_module import Spectrum
from casein_ftir.simulate import simulate_casein, simulate_galalithe


# ---------------------------------------------------------------------------
# Reusable spectra
# ---------------------------------------------------------------------------

@pytest.fixture
def clean_casein_spectrum():
    """Noise-free simulated bovine casein spectrum."""
    return simulate_casein(wn_min=600, wn_max=4000, resolution=2.0,
                            species="bovine", noise_level=0.0)


@pytest.fixture
def noisy_casein_spectrum():
    """Slightly noisy bovine casein spectrum."""
    return simulate_casein(wn_min=600, wn_max=4000, resolution=2.0,
                            species="bovine", noise_level=0.005)


@pytest.fixture
def clean_galalithe_spectrum():
    """Noise-free simulated galalithe spectrum (70% crosslink)."""
    return simulate_galalithe(wn_min=600, wn_max=4000, resolution=2.0,
                                crosslink_degree=0.7, residual_formol=0.02,
                                noise_level=0.0)


@pytest.fixture
def trivial_spectrum():
    """Tiny synthetic spectrum, useful for fast unit checks."""
    wn = np.linspace(4000, 600, 100)  # descending
    ab = np.exp(-((wn - 1648.0) / 30.0) ** 2)  # a single Gaussian
    return Spectrum(wavenumber=wn, absorbance=ab, name="trivial")


@pytest.fixture
def temp_csv(tmp_path):
    """Path to a temp CSV file with a known two-column spectrum."""
    path = tmp_path / "spec.csv"
    wn = np.linspace(4000, 600, 200)
    ab = 0.1 + 0.9 * np.exp(-((wn - 1648.0) / 30.0) ** 2)
    with open(path, "w") as f:
        f.write("# header line 1\n")
        f.write("wavenumber,absorbance\n")
        for w, a in zip(wn, ab):
            f.write(f"{w:.4f},{a:.6f}\n")
    return path
