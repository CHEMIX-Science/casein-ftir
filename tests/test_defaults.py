"""The public package must not silently bundle or fabricate an internal reference."""
from pathlib import Path
import pytest
from casein_ftir import defaults


def test_public_catalog_is_empty():
    assert defaults.list_dft_references() == {}


def test_legacy_galalithe_request_explains_missing_data():
    with pytest.raises(ValueError, match="not included in the public distribution"):
        defaults.get_galalithe_dft_reference()


def test_unknown_reference_is_not_substituted():
    with pytest.raises(KeyError):
        defaults.get_dft_reference("unknown")


def test_internal_reference_file_is_absent():
    assert not (Path(defaults.__file__).parent / "data" / "galalithe_dft_reference.csv").exists()
