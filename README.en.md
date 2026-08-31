# CHEMIX · Casein FTIR

[Documentation française](README.md) · [Scientific scope](docs/SCIENTIFIC_SCOPE.md)

**CONFIRMÉ / confirmed** — Python software developed within CHEMIX for its
Galalithe project. The public preparation includes spectral import, preprocessing,
band integration, reference comparison, Gaussian Amide I fitting, synthetic
examples and local text/HTML reports.

**PROPOSITION / exploratory model** — This is a research and teaching tool.
**À VALIDER / requires validation** — Band assignments, thresholds, structural
interpretations are not independently validated here. The internal DFT
reference is excluded; users can import their own authorized reference. Similarity does not establish purity, safety, crosslink identity,
formaldehyde concentration or suitability for use. Fitted component area fractions
are not experimentally established secondary-structure populations.

[Source code](https://github.com/CHEMIX-Science/casein-ftir) · [Releases](https://github.com/CHEMIX-Science/casein-ftir/releases) · [Automated checks](https://github.com/CHEMIX-Science/casein-ftir/actions)

## Quick start

Python 3.10+; open a terminal in this repository and create a virtual environment:

```bash
python -m venv .venv
```

Activate it on Linux/macOS with `source .venv/bin/activate`. On Windows,
use `.venv\Scripts\python.exe` in place of `python` in the following commands.

```bash
python -m pip install "."
python -m casein_ftir --help
python -m casein_ftir analyze examples/synthetic_casein.csv --deconv --report report.html
python -m casein_ftir compare examples/synthetic_casein.csv examples/synthetic_galalithe.csv
```

**CONFIRMÉ** — Example CSVs are synthetic, not experimental CHEMIX measurements.
CSV input must contain two finite numeric columns: wavenumber (cm⁻¹), absorbance.
Decimal commas and transmittance CSVs are not supported. Optional readers:
`python -m pip install ".[jcamp,opus]"`. Thermo `.spa` is not implemented.

**CONFIRMÉ** — Calibration exists in the Python API, not as a `calibrate` CLI
command. Gaussian conversion reads local output files; it does not run Gaussian.
See [the French guide](README.md) for full commands and limitations.

## Development

```bash
python -m pip install -e ".[test]"
python -m pytest
python -m pip install build
python -m build
```

**CONFIRMÉ / confirmed** — Distributed under the [MIT License](LICENSE), approved
by CHEMIX on 31 August 2026. Collective attribution: CHEMIX contributors.
See [LICENSING.md](LICENSING.md) for scope; scientific validation remains separate.
