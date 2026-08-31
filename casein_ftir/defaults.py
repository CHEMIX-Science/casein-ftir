"""Compatibility accessors for references excluded from the public distribution.

No internal computed or experimental reference is distributed. Users may load
an authorized reference via load_reference() or the CLI --reference option.
"""


def list_dft_references() -> dict:
    """The public catalog is empty; no internal DFT dataset is bundled."""
    return {}


def get_dft_reference(name: str):
    """Keep legacy API calls explicit instead of substituting synthetic data."""
    if name == "galalithe":
        raise ValueError(
            "The internal DFT reference is not included in the public distribution. "
            "Provide your own authorized CSV with --reference or load_reference()."
        )
    raise KeyError(f"Unknown DFT reference {name!r}; the public catalog is empty")


def get_galalithe_dft_reference():
    """Explain the missing internal dataset to callers of the legacy API."""
    return get_dft_reference("galalithe")
