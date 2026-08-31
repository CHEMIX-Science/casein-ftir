"""
Tests for casein_ftir.database

Vérifie l'intégrité interne de la base de données de bandes IR :
- toutes les bandes ont des fenêtres valides englobant leur centre
- toutes les références citées dans les bandes existent dans REFERENCES
- les intensités relatives sont positives
- amide_I a bien rel_intensity_ref = 1.00 (référence)
"""

import pytest

from casein_ftir.database import (
    CASEIN_BANDS,
    AMIDE_I_SUBBANDS,
    GALALITHE_MARKERS,
    REFERENCES,
    IRBand,
    get_band,
    get_reference,
)


ALL_BAND_DICTS = [CASEIN_BANDS, AMIDE_I_SUBBANDS, GALALITHE_MARKERS]


# ---------------------------------------------------------------------------
# Window integrity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("band_dict", ALL_BAND_DICTS)
def test_band_window_contains_center(band_dict):
    """Le centre de chaque bande doit tomber dans sa fenêtre."""
    for name, band in band_dict.items():
        wmin, wmax = band.window
        # exception : free_aldehyde_CO peut avoir un centre en bord de fenêtre
        assert wmin <= band.center <= wmax, (
            f"Band {name}: center {band.center} outside window {band.window}"
        )


@pytest.mark.parametrize("band_dict", ALL_BAND_DICTS)
def test_band_window_well_ordered(band_dict):
    for name, band in band_dict.items():
        wmin, wmax = band.window
        assert wmin < wmax, f"Band {name} has window {band.window}"


@pytest.mark.parametrize("band_dict", ALL_BAND_DICTS)
def test_band_rel_intensity_non_negative(band_dict):
    for name, band in band_dict.items():
        assert band.rel_intensity_ref >= 0.0, (
            f"Band {name} has negative rel_intensity"
        )


# ---------------------------------------------------------------------------
# Reference integrity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("band_dict", ALL_BAND_DICTS)
def test_band_sources_all_resolvable(band_dict):
    """Toute source citée doit être déclarée dans REFERENCES."""
    for name, band in band_dict.items():
        for src in band.sources:
            assert src in REFERENCES, (
                f"Band {name} cites unknown reference {src!r}"
            )


def test_at_least_one_source_per_band():
    for band_dict in ALL_BAND_DICTS:
        for name, band in band_dict.items():
            assert len(band.sources) >= 1, f"Band {name} has no source"


def test_references_have_non_empty_strings():
    for key, ref in REFERENCES.items():
        assert isinstance(ref, str) and len(ref) > 30, (
            f"Reference {key} looks too short or invalid"
        )


# ---------------------------------------------------------------------------
# Amide I reference convention
# ---------------------------------------------------------------------------

def test_amide_I_is_reference():
    """Amide I doit avoir une intensité relative de 1.00 (réf)."""
    assert CASEIN_BANDS["amide_I"].rel_intensity_ref == 1.00


def test_amide_I_center_in_canonical_range():
    """Amide I doit être centré dans la zone 1640-1660 cm-1."""
    c = CASEIN_BANDS["amide_I"].center
    assert 1640.0 <= c <= 1660.0


def test_amide_II_center_in_canonical_range():
    c = CASEIN_BANDS["amide_II"].center
    assert 1510.0 <= c <= 1560.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def test_get_band_returns_correct_object():
    band = get_band("amide_I")
    assert isinstance(band, IRBand)
    assert band.name == "amide_I"


def test_get_band_unknown_raises():
    with pytest.raises(KeyError):
        get_band("not_a_real_band")


def test_get_reference_returns_string():
    ref = get_reference("barth2007")
    assert "Barth" in ref


# ---------------------------------------------------------------------------
# Sub-bands of amide I
# ---------------------------------------------------------------------------

def test_amide_I_subbands_inside_amide_I_window():
    """Toutes les sous-bandes doivent rentrer dans la fenêtre amide I."""
    amide_window = CASEIN_BANDS["amide_I"].window
    for name, sub in AMIDE_I_SUBBANDS.items():
        wmin, wmax = sub.window
        assert wmin >= amide_window[0] and wmax <= amide_window[1], (
            f"Sub-band {name} window {sub.window} not in {amide_window}"
        )


def test_amide_I_subbands_ordered_by_center():
    """Les sous-bandes sont attendues dans un ordre croissant pour faciliter
    la déconvolution."""
    centers = [b.center for b in AMIDE_I_SUBBANDS.values()]
    assert centers == sorted(centers), (
        f"Sub-band centers should be sorted: {centers}"
    )
