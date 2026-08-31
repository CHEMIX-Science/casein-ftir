"""Load the CHEMIX DFT reference shared under the repository MIT license.

The numerical spectrum is inherited unchanged from the internal program.
Its calculation metadata and chemical interpretation remain unverified;
see docs/DATA_PROVENANCE.md. It is not a certified material standard.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict

from .io_module import Spectrum, load_spectrum


# Répertoire data/ relatif au package
_DATA_DIR = Path(__file__).parent / "data"


# Liste des spectres DFT empaquetés
_DFT_REFERENCES: Dict[str, Dict] = {
    "galalithe": {
        "file": "galalithe_dft_reference.csv",
        "description": (
            "Exploratory DFT fragment reference; not a certified material standard. "
            "Legacy calculation metadata (not independently verified): "
            "Model: CH3-CH2-NH-CH2-NH-CH2-CH3 (N,N'-diethylmethanediamine), "
            "minimal model of casein-formaldehyde aminal bridge. "
            "B3LYP/6-311+G(2d,p) + IEFPCM(water), 57 modes, "
            "Legacy scale factor 0.9679 applied; exact provenance À VALIDER. "
            "Bundled with the original CHEMIX program; calculation not independently reproduced."
        ),
        "level": "B3LYP/6-311+G(2d,p)",
        "scale_factor": 0.9679,
        "solvent": "Water (IEFPCM)",
        "n_modes": 57,
        "atoms": 21,
        "molecule": "N,N'-diethylmethanediamine (C5H14N2)",
    },
    # Future : 'casein_phosphate', 'nma_amide_ref', etc.
}


def list_dft_references() -> Dict[str, Dict]:
    """Return the catalog of available DFT reference spectra."""
    return {name: info.copy() for name, info in _DFT_REFERENCES.items()}


def get_dft_reference(name: str) -> Spectrum:
    """Load a packaged DFT reference spectrum.

    Parameters
    ----------
    name : str
        Key from `list_dft_references()`. E.g., 'galalithe'.

    Returns
    -------
    Spectrum with metadata.is_reference = True and
    metadata.type = 'dft_packaged_reference'

    Raises
    ------
    KeyError if `name` is unknown.
    FileNotFoundError if the data file is missing (package install issue).
    """
    if name not in _DFT_REFERENCES:
        raise KeyError(
            f"Unknown DFT reference {name!r}. Available: "
            f"{list(_DFT_REFERENCES.keys())}"
        )
    info = _DFT_REFERENCES[name]
    path = _DATA_DIR / info["file"]
    if not path.exists():
        raise FileNotFoundError(
            f"DFT reference file missing: {path}. "
            "Reinstall the package or recompute via Gaussian."
        )
    spec = load_spectrum(path)
    spec.metadata.update({
        "is_reference":  True,
        "type":          "dft_packaged_reference",
        "reference_id":  name,
        "description":   info["description"],
        "level":         info["level"],
        "scale_factor":  info["scale_factor"],
        "molecule":      info["molecule"],
        "license":       "MIT",
        "scientific_validation": "pending",
    })
    spec.name = f"dft_{name}"
    return spec


def get_galalithe_dft_reference() -> Spectrum:
    """Convenience accessor for the galalithe reference spectrum."""
    return get_dft_reference("galalithe")
