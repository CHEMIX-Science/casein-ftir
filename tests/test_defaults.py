"""Protect the authorized reference values, metadata and fragment interpretation."""
import hashlib
from pathlib import Path
import numpy as np
import pytest
from casein_ftir import defaults
from casein_ftir.reference_module import compare_to_reference


def test_packaged_reference_has_original_values():
    ref = defaults.get_galalithe_dft_reference()
    assert len(ref.wavenumber) == 1801
    np.testing.assert_array_equal(ref.wavenumber, np.arange(4000.0, 399.0, -2.0))
    assert np.isfinite(ref.absorbance).all()
    assert ref.absorbance.min() >= 0.0
    assert ref.absorbance.max() == pytest.approx(1.0)
    assert ref.metadata["is_reference"] is True
    assert ref.metadata["type"] == "dft_packaged_reference"
    assert ref.metadata["license"] == "MIT"
    assert ref.metadata["scientific_validation"] == "pending"


def test_reference_bytes_are_preserved():
    data = (Path(defaults.__file__).parent / "data" / "galalithe_dft_reference.csv").read_bytes()
    assert hashlib.sha256(data).hexdigest() == "ab80a26a2b2b3db9b84ebfe5d4e082b5cf5e117d46fd89323596c9d02d79d3ab"


def test_catalog_cannot_be_mutated_by_callers():
    catalog = defaults.list_dft_references()
    assert set(catalog) == {"galalithe"}
    catalog["galalithe"]["file"] = "other.csv"
    assert defaults.list_dft_references()["galalithe"]["file"] == "galalithe_dft_reference.csv"


def test_unknown_reference_is_not_substituted():
    with pytest.raises(KeyError):
        defaults.get_dft_reference("unknown")


def test_missing_reference_reports_installation_problem(monkeypatch, tmp_path):
    monkeypatch.setattr(defaults, "_DATA_DIR", tmp_path)
    with pytest.raises(FileNotFoundError, match="Reinstall"):
        defaults.get_galalithe_dft_reference()


def test_reference_comparison_does_not_grade_material_quality():
    ref = defaults.get_galalithe_dft_reference()
    comparison = compare_to_reference(ref, ref)
    assert comparison["overall_quality"] == "fragment_reference"
    assert "do not establish" in comparison["recommended_action"]
