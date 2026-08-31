"""
casein_ftir.io_module
=====================

Lecture de spectres IR depuis différents formats :
    - CSV / TSV / TXT (deux colonnes : wavenumber, absorbance)
    - JCAMP-DX (.dx, .jdx, .jcm)
    - Bruker OPUS (binaire, optionnel)
    - Thermo .spa non pris en charge : exporter en CSV

Tous les readers retournent un objet `Spectrum` standardisé.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, Any
import numpy as np
import csv


# ---------------------------------------------------------------------------
# Standard data structure
# ---------------------------------------------------------------------------

@dataclass
class Spectrum:
    """Standardized FTIR spectrum container.

    Convention : wavenumber sorted in DESCENDING order (common IR usage:
    4000 -> 400 cm^-1). Absorbance is unitless (or in absorbance units, AU).
    """
    wavenumber: np.ndarray              # cm^-1, descending
    absorbance: np.ndarray              # same shape
    name: str = "spectrum"
    source_file: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.wavenumber = np.asarray(self.wavenumber, dtype=float)
        self.absorbance = np.asarray(self.absorbance, dtype=float)
        if self.wavenumber.shape != self.absorbance.shape:
            raise ValueError(
                f"Shape mismatch: wavenumber {self.wavenumber.shape}, "
                f"absorbance {self.absorbance.shape}"
            )
        if self.wavenumber.ndim != 1:
            raise ValueError("Spectrum arrays must be one-dimensional")
        if not (np.isfinite(self.wavenumber).all() and np.isfinite(self.absorbance).all()):
            raise ValueError("Spectrum values must be finite (no NaN or infinity)")
        order = np.argsort(self.wavenumber)[::-1]
        self.wavenumber = self.wavenumber[order]
        self.absorbance = self.absorbance[order]
        if len(self.wavenumber) > 1 and np.any(np.diff(self.wavenumber) == 0):
            raise ValueError("Duplicate wavenumbers: aggregate duplicates explicitly before import")

    @property
    def n_points(self) -> int:
        return len(self.wavenumber)

    @property
    def range(self) -> tuple:
        return float(self.wavenumber.min()), float(self.wavenumber.max())

    def slice(self, wn_min: float, wn_max: float) -> "Spectrum":
        """Return a sub-spectrum restricted to [wn_min, wn_max]."""
        mask = (self.wavenumber >= wn_min) & (self.wavenumber <= wn_max)
        return Spectrum(
            wavenumber=self.wavenumber[mask],
            absorbance=self.absorbance[mask],
            name=f"{self.name}[{wn_min}-{wn_max}]",
            source_file=self.source_file,
            metadata=dict(self.metadata),
        )

    def copy(self) -> "Spectrum":
        return Spectrum(
            wavenumber=self.wavenumber.copy(),
            absorbance=self.absorbance.copy(),
            name=self.name,
            source_file=self.source_file,
            metadata=dict(self.metadata),
        )

    def __repr__(self) -> str:
        wn0, wn1 = self.range
        return (
            f"Spectrum(name={self.name!r}, n={self.n_points}, "
            f"range={wn1:.1f}-{wn0:.1f} cm-1)"
        )


# ---------------------------------------------------------------------------
# Format-specific readers
# ---------------------------------------------------------------------------

def _read_csv(path: Path) -> Spectrum:
    """Read exactly two numeric columns; never silently discard malformed data."""
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    metadata = {}
    for line in lines:
        if line.startswith("# type: "):
            kind = line.partition(": ")[2].strip()
            if kind in {"simulated", "dft_simulated", "dft_packaged_reference", "difference", "subtracted"}:
                metadata["type"] = kind
    content = [line for line in lines if line.strip() and not line.lstrip().startswith("#")]
    for delimiter in (",", "\t", ";", None):
        rows = []
        headers = 0
        failed = False
        for line in content:
            fields = line.split() if delimiter is None else next(csv.reader([line], delimiter=delimiter))
            try:
                if len(fields) != 2:
                    raise ValueError("Expected exactly two columns")
                row = [float(value) for value in fields]
            except ValueError:
                # Accept up to five plain-text instrument/header rows before data.
                numeric_field = False
                for value in fields:
                    try:
                        float(value)
                        numeric_field = True
                    except ValueError:
                        pass
                if rows or headers >= 5 or numeric_field:
                    failed = True
                    break
                headers += 1
                continue
            rows.append(row)
        if not failed and len(rows) >= 2:
            arr = np.asarray(rows)
            return Spectrum(arr[:, 0], arr[:, 1], name=path.stem,
                            source_file=str(path), metadata=metadata)
    raise ValueError(f"Cannot parse {path.name}: use two finite numeric columns, "
                     "a decimal point, and comma/tab/semicolon/whitespace separators")


def _read_jcamp(path: Path) -> Spectrum:
    """Read a JCAMP-DX file using the `jcamp` package."""
    try:
        from jcamp import jcamp_readfile
    except ImportError as e:
        raise ImportError(
            "Reading JCAMP-DX requires the 'jcamp' package. "
            "Install with: pip install jcamp"
        ) from e
    data = jcamp_readfile(str(path))
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    # If file was in transmittance, convert to absorbance:
    yunits = (data.get("yunits") or "").upper()
    if "TRANSMITTANCE" in yunits:
        # protect against zero/negative T
        y = np.clip(y, 1e-6, None)
        # Normalise transmittance to [0,1] if expressed in percent
        if y.max() > 1.5:
            y = y / 100.0
        y = -np.log10(y)
    return Spectrum(
        wavenumber=x,
        absorbance=y,
        name=path.stem,
        source_file=str(path),
        metadata={k: data.get(k) for k in ("title", "owner", "origin", "date")},
    )


def _read_opus(path: Path) -> Spectrum:
    """Read a Bruker OPUS binary file (optional dependency)."""
    try:
        from brukeropusreader import read_file
    except ImportError as e:
        raise ImportError(
            "Reading Bruker OPUS requires the 'brukeropusreader' package. "
            "Install with: pip install brukeropusreader"
        ) from e
    opus = read_file(str(path))
    # Pick the absorbance block in priority, then transmittance
    for key in ("AB", "Absorbance"):
        if key in opus:
            ab = np.asarray(opus[key], dtype=float)
            x = opus.get_range(key)
            return Spectrum(
                wavenumber=np.asarray(x, dtype=float),
                absorbance=ab,
                name=path.stem,
                source_file=str(path),
                metadata={"format": "OPUS"},
            )
    for key in ("TR", "Transmittance"):
        if key in opus:
            tr = np.clip(np.asarray(opus[key], dtype=float), 1e-6, None)
            if tr.max() > 1.5:
                tr = tr / 100.0
            ab = -np.log10(tr)
            x = opus.get_range(key)
            return Spectrum(
                wavenumber=np.asarray(x, dtype=float),
                absorbance=ab,
                name=path.stem,
                source_file=str(path),
                metadata={"format": "OPUS", "converted_from": "transmittance"},
            )
    raise IOError(f"No AB or TR block found in OPUS file: {path}")


# ---------------------------------------------------------------------------
# Public dispatcher with format auto-detection
# ---------------------------------------------------------------------------

_EXT_DISPATCH = {
    ".csv": _read_csv,
    ".tsv": _read_csv,
    ".txt": _read_csv,
    ".dat": _read_csv,
    ".asc": _read_csv,
    ".dx":  _read_jcamp,
    ".jdx": _read_jcamp,
    ".jcm": _read_jcamp,
    ".dpt": _read_csv,  # Bruker plain text export
}


def load_spectrum(path: str | Path, fmt: Optional[str] = None) -> Spectrum:
    """Load a spectrum from disk, auto-detecting format from extension.

    Parameters
    ----------
    path : str or Path
        Path to the spectrum file.
    fmt : str, optional
        Force a format among {'csv', 'jcamp', 'opus'}. Default: auto-detect.

    Returns
    -------
    Spectrum
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if fmt is not None:
        fmt = fmt.lower()
        if fmt == "csv":
            return _read_csv(path)
        if fmt == "jcamp":
            return _read_jcamp(path)
        if fmt == "opus":
            return _read_opus(path)
        raise ValueError(f"Unknown forced format: {fmt!r}")

    ext = path.suffix.lower()
    if ext == ".spa":
        raise ValueError("Thermo .spa is not supported; export wavenumber/absorbance to CSV")
    if ext in _EXT_DISPATCH:
        return _EXT_DISPATCH[ext](path)

    # Files without standard extension : try OPUS (binary), then CSV
    try:
        return _read_opus(path)
    except Exception:
        return _read_csv(path)


def save_spectrum_csv(spec: Spectrum, path: str | Path) -> None:
    """Export a Spectrum to a simple CSV file."""
    path = Path(path)
    with path.open("w", encoding="utf-8") as f:
        f.write("# wavenumber_cm-1,absorbance\n")
        f.write(f"# name: {spec.name}\n")
        if spec.source_file:
            f.write(f"# source: {Path(spec.source_file).name}\n")
        if spec.metadata.get("type"):
            f.write(f"# type: {spec.metadata['type']}\n")
        for wn, a in zip(spec.wavenumber, spec.absorbance):
            f.write(f"{wn:.4f},{a:.6f}\n")
