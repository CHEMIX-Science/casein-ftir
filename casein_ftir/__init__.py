"""
casein_ftir
===========

Analyse quantitative de spectres IR pour l'étude caséine -> galalithe.

Modules principaux :
    - io_module     : lecture multi-format
    - preprocess    : baseline, smoothing, normalisation
    - database      : positions de bandes IR (litterature)
    - peak_detection: detection automatique
    - quantification: aires, ratios, Beer-Lambert
    - amide_deconv  : deconvolution amide I
    - compare       : comparaison avant/apres
    - simulate      : spectres theoriques
    - report        : export rapports
"""

__version__ = "0.2.0"

from .io_module import Spectrum, load_spectrum
from .database import CASEIN_BANDS, AMIDE_I_SUBBANDS, GALALITHE_MARKERS
from .reference_module import (load_reference, compare_to_reference,
                                  subtract_reference)
from .defaults import (get_dft_reference, get_galalithe_dft_reference,
                          list_dft_references)

__all__ = [
    "Spectrum",
    "load_spectrum",
    "CASEIN_BANDS",
    "AMIDE_I_SUBBANDS",
    "GALALITHE_MARKERS",
    "load_reference",
    "compare_to_reference",
    "subtract_reference",
    "get_dft_reference",
    "get_galalithe_dft_reference",
    "list_dft_references",
]
