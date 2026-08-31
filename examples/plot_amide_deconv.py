"""
plot_amide_deconv.py
====================

Standalone helper script to plot the Amide I deconvolution from a CSV
exported by :

    casein-ftir deconv input.csv --export deconv.csv

Usage :
    python plot_amide_deconv.py deconv.csv
    python plot_amide_deconv.py deconv.csv --output figure.png
    python plot_amide_deconv.py deconv.csv --output figure.pdf --dpi 300

Generates a publication-style figure with :
    - the baseline-corrected observed Amide I band (points)
    - the total Gaussian fit (line, dashed)
    - the 6 individual sub-bands, color-coded with labels and %
    - the residual at the bottom

Requires : numpy, pandas, matplotlib (installed with the package).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Couleurs assignées à chaque sous-bande (cohérentes avec la littérature)
SUBBAND_COLORS = {
    "beta_aggregate":    "#8b0000",   # rouge foncé
    "beta_sheet":        "#cc4444",   # rouge
    "random_coil":       "#888888",   # gris
    "alpha_helix":       "#4477aa",   # bleu
    "beta_turn":         "#aa6633",   # marron
    "beta_antiparallel": "#dd7788",   # rose
}

SUBBAND_LABELS = {
    "beta_aggregate":    "β-aggregate",
    "beta_sheet":        "β-sheet",
    "random_coil":       "random coil",
    "alpha_helix":       "α-helix",
    "beta_turn":         "β-turn",
    "beta_antiparallel": "β-antiparallel",
}


def parse_metadata(path: Path) -> dict:
    """Read the % comments at the top of the file to retrieve sub-band %."""
    info = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break
            line = line.strip()[1:].strip()
            if "pct=" in line:
                # Example : 'random_coil          center=1644.0 ... pct=36.21%'
                parts = line.split()
                name = parts[0]
                pct_str = [p for p in parts if p.startswith("pct=")][0]
                pct = float(pct_str.replace("pct=", "").replace("%", ""))
                info[name] = pct
            elif "R^2" in line:
                info["r2"] = float(line.split("=")[1].strip())
    return info


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", help="CSV file produced by casein-ftir --deconv-export")
    p.add_argument("--output", "-o", help="output figure (.png, .pdf, .svg)")
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--show", action="store_true",
                   help="show figure interactively (default if no --output)")
    args = p.parse_args()

    path = Path(args.csv)
    if not path.exists():
        sys.exit(f"File not found: {path}")

    df = pd.read_csv(path, comment="#")
    meta = parse_metadata(path)

    x = df["wavenumber"].values
    obs = df["observed"].values
    fit = df["fit_total"].values
    res = df["residual"].values

    # Layout : main plot + residual below
    fig, (ax, ax_res) = plt.subplots(
        2, 1, figsize=(8, 6), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05},
    )

    # ---- Main plot ----
    ax.scatter(x, obs, s=12, color="black", alpha=0.6, label="Observed",
                zorder=3)
    ax.plot(x, fit, color="black", linestyle="--", linewidth=1.5,
             label=f"Total fit (R²={meta.get('r2', 0):.4f})", zorder=4)

    sub_band_names = ["beta_aggregate", "beta_sheet", "random_coil",
                       "alpha_helix", "beta_turn", "beta_antiparallel"]
    for name in sub_band_names:
        if name not in df.columns:
            continue
        pct = meta.get(name, 0.0)
        ax.fill_between(x, 0, df[name].values,
                          color=SUBBAND_COLORS[name], alpha=0.3,
                          label=f"{SUBBAND_LABELS[name]} ({pct:.1f}%)")
        ax.plot(x, df[name].values, color=SUBBAND_COLORS[name],
                 linewidth=1.2, alpha=0.85)

    ax.set_ylabel("Absorbance (a.u.)", fontsize=11)
    ax.set_title("CHEMIX — exploratory Amide I fit\nModel area fractions, not validated structural populations",
                  fontsize=12)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.invert_xaxis()   # IR convention : descending wavenumber
    ax.grid(alpha=0.3)

    # ---- Residual plot ----
    ax_res.plot(x, res, color="black", linewidth=0.8)
    ax_res.axhline(0, color="gray", linewidth=0.5)
    ax_res.set_xlabel("Wavenumber (cm⁻¹)", fontsize=11)
    ax_res.set_ylabel("Residual", fontsize=10)
    ax_res.grid(alpha=0.3)


    if args.output:
        fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
        print(f"Figure saved to {args.output}")
    if args.show or not args.output:
        plt.show()


if __name__ == "__main__":
    main()
