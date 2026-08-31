"""
casein_ftir.database
====================

Base de données des bandes IR caractéristiques de la caséine.

Chaque entrée fournit :
    - center      : position centrale (cm^-1)
    - window      : intervalle [min, max] dans lequel chercher le pic (cm^-1)
    - assignment  : attribution vibrationnelle
    - rel_intensity_ref : intensité relative typique (amide I = 1.00)
    - source      : référence bibliographique (clé du dict REFERENCES)

PROPOSITION — Fenêtres, intensités et attributions héritées du modèle interne.
À VALIDER — Traçabilité de chaque valeur et validité sur chaque échantillon.
Les références bibliographiques ne certifient pas chaque paramètre numérique.
"""

from dataclasses import dataclass, field
from typing import Tuple, Dict, List


# ---------------------------------------------------------------------------
# Références bibliographiques
# ---------------------------------------------------------------------------

REFERENCES: Dict[str, str] = {
    "barth2007": (
        "Barth A. (2007). Infrared spectroscopy of proteins. "
        "Biochim. Biophys. Acta - Bioenergetics 1767(9), 1073-1101. "
        "DOI: 10.1016/j.bbabio.2007.06.004"
    ),
    "dimitriou2025": (
        "Dimitriou L. et al. (2025). Rapid Classification of Cow, Goat, "
        "and Sheep Milk Using ATR-FTIR and Multivariate Analysis. "
        "Sci 7(3), 87. DOI: 10.3390/sci7030087"
    ),
    "lecaroz2025": (
        "D'Incecco P. et al. (2025). Functionalization of Sodium "
        "Caseinate for Production of Neat Films: Effects of Casein Crosslinking Induced by Heating at Alkaline pH or Light Exposure. Foods 14(16), 2764. "
        "DOI: 10.3390/foods14162764"
    ),
    "demeutter2020": (
        "Sadat A., Joye I.J. (2020). Peak Fitting Applied to "
        "FTIR and Raman Spectroscopic Analysis of Proteins. "
        "Appl. Sci. 10(17), 5918. DOI: 10.3390/app10175918"
    ),
    "cirak2018": (
        "Cirak O., Icyer N.C., Durak M.Z. (2018). Rapid detection of "
        "adulteration of milks from different species using FTIR. "
        "J. Dairy Res. 85(2), 222-225. DOI: 10.1017/S0022029918000201"
    ),
    "lu2020": (
        "Ji Y. et al. (2020). DFT-Calculated IR Spectrum Amide I, II, and III "
        "Band Contributions of N-Methylacetamide Fine Components. "
        "ACS Omega 5, 8572-8578. DOI: 10.1021/acsomega.9b04421"
    ),
}


# ---------------------------------------------------------------------------
# Définition d'une bande IR
# ---------------------------------------------------------------------------

@dataclass
class IRBand:
    name: str
    center: float                  # cm^-1
    window: Tuple[float, float]    # (min, max) for peak search
    assignment: str
    rel_intensity_ref: float       # relative to amide I = 1.00
    sources: List[str] = field(default_factory=list)
    status: str = "À VALIDER"

    def __str__(self) -> str:
        return f"{self.name} @ {self.center} cm-1 ({self.assignment})"


# ---------------------------------------------------------------------------
# Bandes universelles de la caséine (bovine, caprine, ovine -- très proches)
# ---------------------------------------------------------------------------

CASEIN_BANDS: Dict[str, IRBand] = {
    "amide_A": IRBand(
        name="amide_A",
        center=3288.0,
        window=(3270.0, 3500.0),
        assignment="nu(N-H) elongation + nu(O-H)",
        rel_intensity_ref=0.70,
        sources=["barth2007", "lecaroz2025", "lu2020"],
    ),
    "CH3_asym": IRBand(
        name="CH3_asym",
        center=2960.0,
        window=(2950.0, 2975.0),
        assignment="nu(C-H) asym CH3 aliphatique",
        rel_intensity_ref=0.20,
        sources=["lecaroz2025"],
    ),
    "CH2_asym": IRBand(
        name="CH2_asym",
        center=2926.0,
        window=(2915.0, 2940.0),
        assignment="nu(C-H) asym CH2 aliphatique",
        rel_intensity_ref=0.25,
        sources=["lecaroz2025"],
    ),
    "CH_sym": IRBand(
        name="CH_sym",
        center=2870.0,
        window=(2855.0, 2885.0),
        assignment="nu(C-H) sym CH2/CH3",
        rel_intensity_ref=0.15,
        sources=["lecaroz2025"],
    ),
    "amide_I": IRBand(
        name="amide_I",
        center=1648.0,
        window=(1600.0, 1700.0),
        assignment="nu(C=O) 70-85% + nu(C-N) 10-20% (peptide bond)",
        rel_intensity_ref=1.00,
        sources=["barth2007", "lecaroz2025", "demeutter2020"],
    ),
    "amide_II": IRBand(
        name="amide_II",
        center=1540.0,
        window=(1510.0, 1580.0),
        assignment="delta(N-H) + nu(C-N) peptide bond",
        rel_intensity_ref=0.65,
        sources=["barth2007", "demeutter2020"],
    ),
    "CH2_bend": IRBand(
        name="CH2_bend",
        center=1450.0,
        window=(1430.0, 1470.0),
        assignment="delta(CH2)/delta(CH3) ciseaux",
        rel_intensity_ref=0.30,
        sources=["lecaroz2025"],
    ),
    "COO_sym": IRBand(
        name="COO_sym",
        center=1400.0,
        window=(1385.0, 1415.0),
        assignment="nu sym(COO-) chaines laterales Asp/Glu",
        rel_intensity_ref=0.32,
        sources=["barth2007"],
    ),
    "amide_III": IRBand(
        name="amide_III",
        center=1240.0,
        window=(1230.0, 1300.0),
        assignment="mix nu(C-N) + delta(N-H)",
        rel_intensity_ref=0.20,
        sources=["barth2007", "lu2020"],
    ),
    "phosphate": IRBand(
        name="phosphate",
        center=1100.0,
        window=(1050.0, 1150.0),
        assignment="nu(P=O)/nu(P-O) phosphoserine -- attribution non specifique, a valider",
        rel_intensity_ref=0.40,
        sources=["dimitriou2025"],
    ),
    "C_O_skeletal": IRBand(
        name="C_O_skeletal",
        center=1020.0,
        window=(970.0, 1080.0),
        assignment="nu(C-O) squelette + phosphate (overlap)",
        rel_intensity_ref=0.30,
        sources=["dimitriou2025"],
    ),
}


# ---------------------------------------------------------------------------
# Sous-bandes amide I pour deconvolution -- analyse de structure secondaire
# ---------------------------------------------------------------------------

AMIDE_I_SUBBANDS: Dict[str, IRBand] = {
    "beta_aggregate": IRBand(
        name="beta_aggregate",
        center=1618.0,
        window=(1610.0, 1625.0),
        assignment="feuillets beta intermoleculaires / agregats",
        rel_intensity_ref=0.10,
        sources=["demeutter2020", "barth2007"],
    ),
    "beta_sheet": IRBand(
        name="beta_sheet",
        center=1632.0,
        window=(1625.0, 1640.0),
        assignment="feuillets beta natifs",
        rel_intensity_ref=0.15,
        sources=["demeutter2020", "barth2007"],
    ),
    "random_coil": IRBand(
        name="random_coil",
        center=1644.0,
        window=(1640.0, 1648.0),
        assignment="random coil / desordonne (dominant dans caseine)",
        rel_intensity_ref=0.40,
        sources=["demeutter2020", "barth2007"],
    ),
    "alpha_helix": IRBand(
        name="alpha_helix",
        center=1654.0,
        window=(1648.0, 1660.0),
        assignment="helices alpha",
        rel_intensity_ref=0.20,
        sources=["demeutter2020", "barth2007"],
    ),
    "beta_turn": IRBand(
        name="beta_turn",
        center=1670.0,
        window=(1660.0, 1680.0),
        assignment="tours beta",
        rel_intensity_ref=0.10,
        sources=["demeutter2020", "barth2007"],
    ),
    "beta_antiparallel": IRBand(
        name="beta_antiparallel",
        center=1685.0,
        window=(1680.0, 1695.0),
        assignment="feuillets beta antiparalleles",
        rel_intensity_ref=0.05,
        sources=["demeutter2020", "barth2007"],
    ),
}


# ---------------------------------------------------------------------------
# Bandes additionnelles pour la GALALITHE (caseine + formaldehyde)
# Apparaissent ou s'intensifient lors de la reticulation
# ---------------------------------------------------------------------------

GALALITHE_MARKERS: Dict[str, IRBand] = {
    "methylene_bridge": IRBand(
        name="methylene_bridge",
        center=2920.0,
        window=(2900.0, 2940.0),
        assignment="region C-H exploratoire; attribution aux ponts a valider",
        rel_intensity_ref=0.40,
        sources=["barth2007"],
    ),
    "free_aldehyde_CO": IRBand(
        name="free_aldehyde_CO",
        center=1720.0,
        window=(1700.0, 1740.0),
        assignment="region carbonyle exploratoire; non specifique du formaldehyde",
        rel_intensity_ref=0.0,
        sources=["barth2007"],
    ),
    "N_CH2_N": IRBand(
        name="N_CH2_N",
        center=1370.0,
        window=(1350.0, 1390.0),
        assignment="region exploratoire; attribution N-CH2-N a valider",
        rel_intensity_ref=0.20,
        sources=["barth2007"],
    ),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_band(name: str) -> IRBand:
    """Get an IRBand by name from any of the databases."""
    for db in (CASEIN_BANDS, AMIDE_I_SUBBANDS, GALALITHE_MARKERS):
        if name in db:
            return db[name]
    raise KeyError(f"Unknown band: {name!r}")


def get_reference(key: str) -> str:
    """Return the full bibliographic reference for a source key."""
    return REFERENCES[key]


def print_all_bands() -> None:
    """Print all known bands with their sources (for CLI 'info' command)."""
    print("=== CASEIN BANDS ===")
    for b in CASEIN_BANDS.values():
        srcs = ", ".join(b.sources)
        print(f"  {b.name:20s} {b.center:6.1f} cm-1  [{srcs}] — {b.status}")
    print("\n=== AMIDE I SUB-BANDS ===")
    for b in AMIDE_I_SUBBANDS.values():
        srcs = ", ".join(b.sources)
        print(f"  {b.name:20s} {b.center:6.1f} cm-1  [{srcs}] — {b.status}")
    print("\n=== GALALITHE MARKERS ===")
    for b in GALALITHE_MARKERS.values():
        srcs = ", ".join(b.sources)
        print(f"  {b.name:20s} {b.center:6.1f} cm-1  [{srcs}] — {b.status}")
    print("\n=== REFERENCES ===")
    for key, ref in REFERENCES.items():
        print(f"  [{key}] {ref}")
