# CHEMIX · Casein FTIR

**A Python toolkit for exploring infrared spectra of casein.**

Developed by CHEMIX for the Galalithe project, Casein FTIR takes you from spectral
import to band comparison and a local analysis report. It is intended for
students, educators and researchers working with casein and related materials.

[Français](README.md) · [Downloads](https://github.com/CHEMIX-Science/casein-ftir/releases)

## What it does

- Correct baselines, smooth signals and normalize spectra.
- Find peaks in predefined band windows, integrate areas and calculate ratios.
- Fit six Gaussian components in the Amide I region and plot the decomposition.
- Compare spectra, examine band differences and export numerical differences.
- Use the included Galalithe DFT fragment reference or load your own reference.
- Convert Gaussian frequency-calculation output into an IR spectrum.
- Save text or HTML reports locally and try reproducible synthetic examples.

This is an exploratory tool. Similarity scores and fitted components alone do
not establish material identity, purity, safety or crosslinking degree. The DFT
reference models a small molecular fragment, not an entire protein or a certified
material standard. See [methods and limitations](docs/SCIENTIFIC_SCOPE.md).

## Quick start

Install Python 3.10+, download and unpack the source, then run from its folder:

```bash
python -m venv .venv
```

Activate it on Linux/macOS with `source .venv/bin/activate`. On Windows, use
`.venv\Scripts\python.exe` in place of `python` below.

```bash
python -m pip install "."
python -m casein_ftir analyze examples/synthetic_casein.csv --deconv --deconv-export deconv.csv --report report.html
python examples/plot_amide_deconv.py deconv.csv --output figure.png
```

Open `report.html` and `figure.png` to inspect the result. The example CSVs are
synthetic, not experimental measurements. Replace the input path with your own
spectrum and choose distinct output names to avoid overwriting files.

Compare against the included theoretical reference:

```bash
python -m casein_ftir analyze examples/synthetic_galalithe.csv --galalithe --use-default-galalithe-ref --no-preprocess-ref --report dft-comparison.html
```

For your own reference, use `--reference my_reference.csv`. Run
`python -m casein_ftir info` to see the bundled reference description.

## Input and extensions

CSV input uses two finite numeric columns: wavenumber in cm⁻¹ and absorbance.
Use decimal points; convert transmittance to absorbance before importing.
JCAMP-DX and Bruker OPUS readers are optional (`pip install ".[jcamp,opus]"`);
check compatibility with your instrument files. Thermo `.spa` is not supported.

Optional fitting and baseline engines: `pip install ".[deconv,baseline]"`.
Gaussian conversion reads calculation output; it does not run Gaussian.
See [the French guide](README.md) for commands, formats and examples.

## Development and license

```bash
python -m pip install -e ".[test]"
python -m pytest
```

[Testing](docs/VALIDATION.md) · [Contributing](CONTRIBUTING.md) · [Data](docs/DATA_PROVENANCE.md)

The code, documentation, examples and bundled DFT reference are provided under
the [MIT License](LICENSE). See [attribution and scope](LICENSING.md).

Developed by **[CHEMIX](https://chemix-paris.com)** for **Galalithe**.
