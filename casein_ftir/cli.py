"""
casein_ftir.cli
===============

Interface en ligne de commande.

Usage :
    casein-ftir <command> [options]

Commandes :
    analyze    : analyse complète d'un spectre
    compare    : comparaison avant/après de deux spectres
    deconv     : déconvolution amide I uniquement
    simulate   : génère un spectre théorique
    info       : liste les bandes connues et les références
"""

from __future__ import annotations
import argparse
import sys
import json
from pathlib import Path

from . import __version__
from .io_module import load_spectrum, save_spectrum_csv
from .preprocess import full_preprocess
from .peak_detection import detect_bands, report_peaks
from .quantification import diagnostic_ratios, galalithe_indicators
from .amide_deconv import deconvolve_amide_I, report_deconvolution
from .compare import (difference_spectrum, highlight_changes,
                       similarity_metric, report_changes)
from .simulate import simulate_casein, simulate_galalithe
from .database import print_all_bands
from .report import text_report, html_report


# ---------------------------------------------------------------------------
# Sub-commands
# ---------------------------------------------------------------------------

def cmd_analyze(args):
    spec = load_spectrum(args.input, fmt=args.format)
    print(f"Loaded: {spec!r}")

    if not args.no_preprocess:
        spec = full_preprocess(
            spec,
            baseline=args.baseline,
            smooth=not args.no_smooth,
            normalize=args.normalize,
        )
        print(f"Preprocessed: {spec.metadata}")

    peaks = detect_bands(spec, band_set="all")
    ratios = diagnostic_ratios(spec)
    galalithe = galalithe_indicators(spec) if args.galalithe else None

    deconv = None
    if args.deconv:
        deconv = deconvolve_amide_I(spec)
        if args.deconv_export:
            from .amide_deconv import export_deconvolution_csv
            export_deconvolution_csv(deconv, args.deconv_export)
            print(f"Amide I deconvolution exported to {args.deconv_export}")

    # --- Reference comparison (optional) -----------------------------------
    reference_analysis = None
    ref = None
    if args.use_default_galalithe_ref:
        from .defaults import get_galalithe_dft_reference
        ref = get_galalithe_dft_reference()
        print(f"Using packaged DFT galalithe reference: {ref.name}")
        print(f"  Description: {ref.metadata['description'][:80]}...")
    elif args.reference:
        from .reference_module import load_reference
        ref = load_reference(
            args.reference,
            preprocess=not args.no_preprocess_ref,
            baseline=args.baseline,
            smooth=not args.no_smooth,
            normalize=args.normalize,
        )
        print(f"Reference loaded: {ref!r}")

    if ref is not None:
        from .reference_module import compare_to_reference, subtract_reference
        reference_analysis = compare_to_reference(
            spec, ref, band_change_threshold=args.ref_threshold
        )
        if args.subtract_ref:
            diff = subtract_reference(spec, ref, scale=args.subtract_scale)
            from .io_module import save_spectrum_csv
            save_spectrum_csv(diff, args.subtract_ref)
            print(f"Sample - reference saved to {args.subtract_ref}")

    report = text_report(spec, peaks=peaks, ratios=ratios,
                          deconv=deconv, galalithe=galalithe,
                          reference_analysis=reference_analysis)
    print(report)

    if args.report:
        path = Path(args.report)
        if path.suffix.lower() == ".html":
            html_report(spec, peaks=peaks, ratios=ratios,
                         deconv=deconv, galalithe=galalithe,
                         reference_analysis=reference_analysis,
                         output=str(path))
        else:
            text_report(spec, peaks=peaks, ratios=ratios,
                         deconv=deconv, galalithe=galalithe,
                         reference_analysis=reference_analysis,
                         output=str(path))
        print(f"Report saved to {path}")


def cmd_compare(args):
    before = load_spectrum(args.before)
    after  = load_spectrum(args.after)
    print(f"Before: {before!r}")
    print(f"After : {after!r}")

    if not args.no_preprocess:
        before = full_preprocess(before, baseline="als", normalize="amide_I")
        after  = full_preprocess(after,  baseline="als", normalize="amide_I")

    changes = highlight_changes(before, after, threshold=args.threshold)
    sim = similarity_metric(before, after)

    print(report_changes(changes))
    print("\nSimilarity metrics:")
    for k, v in sim.items():
        print(f"  {k:20s} = {v:.4f}")

    if args.diff_out:
        diff = difference_spectrum(before, after, normalize=False)
        save_spectrum_csv(diff, args.diff_out)
        print(f"Difference spectrum saved to {args.diff_out}")


def cmd_deconv(args):
    spec = load_spectrum(args.input)
    if not args.no_preprocess:
        spec = full_preprocess(spec, baseline="als", normalize="amide_I")
    result = deconvolve_amide_I(
        spec, wn_min=args.wn_min, wn_max=args.wn_max
    )
    print(report_deconvolution(result))
    if args.export:
        from .amide_deconv import export_deconvolution_csv
        export_deconvolution_csv(result, args.export)
        print(f"\nDeconvolution data exported to {args.export}")


def cmd_simulate(args):
    if args.kind == "casein":
        spec = simulate_casein(
            wn_min=args.wn_min, wn_max=args.wn_max,
            species=args.species, noise_level=args.noise,
        )
    else:
        spec = simulate_galalithe(
            wn_min=args.wn_min, wn_max=args.wn_max,
            crosslink_degree=args.crosslink,
            residual_formol=args.formol,
            species=args.species, noise_level=args.noise,
        )
    save_spectrum_csv(spec, args.output)
    print(f"Simulated spectrum written to {args.output}")
    print(f"  metadata: {spec.metadata}")


def cmd_info(args):
    print(f"casein_ftir v{__version__}")
    print()
    print_all_bands()
    # also list DFT references
    from .defaults import list_dft_references
    refs = list_dft_references()
    if not refs:
        print("\nNo internal DFT reference is bundled. Use --reference with an authorized CSV.")
    if refs:
        print("\n=== PACKAGED DFT REFERENCES ===")
        for name, info in refs.items():
            print(f"  [{name}]")
            print(f"    molecule : {info['molecule']}")
            print(f"    level    : {info['level']}")
            print(f"    solvent  : {info['solvent']}")
            print(f"    modes    : {info['n_modes']} ({info['atoms']} atoms)")
            print(f"    desc     : {info['description'][:100]}...")


def cmd_gaussian(args):
    from .io_gaussian import (gaussian_log_to_spectrum, list_modes,
                                SCALE_FACTORS)
    if args.action == "list-scales":
        print("Legacy scale-factor presets — À VALIDER against the exact calculation level:")
        for k, v in SCALE_FACTORS.items():
            print(f"  {k:30s} -> {v:.4f}")
        return

    if args.action == "inspect":
        text = list_modes(args.input, scale=args.scale, top_n=args.top_n)
        print(text)
        return

    if args.action == "convert":
        spec = gaussian_log_to_spectrum(
            args.input,
            scale=args.scale,
            level=args.level,
            wn_min=args.wn_min,
            wn_max=args.wn_max,
            resolution=args.resolution,
            fwhm=args.fwhm,
            lineshape=args.lineshape,
        )
        save_spectrum_csv(spec, args.output)
        print(f"Gaussian-derived spectrum written to {args.output}")
        print(f"  metadata: {spec.metadata}")
        return

    raise ValueError(f"Unknown gaussian action: {args.action}")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="casein-ftir",
        description="CHEMIX — analyse FTIR exploratoire de caséine et galalithe",
    )
    p.add_argument("--version", action="version",
                   version=f"casein_ftir {__version__}")

    sub = p.add_subparsers(dest="command", required=True)

    # analyze
    pa = sub.add_parser("analyze", help="full analysis of one spectrum")
    pa.add_argument("input", help="path to spectrum file")
    pa.add_argument("--format", choices=["csv", "jcamp", "opus"], default=None)
    pa.add_argument("--no-preprocess", action="store_true")
    pa.add_argument("--baseline", choices=["als", "poly", "rubberband"],
                    default="als")
    pa.add_argument("--no-smooth", action="store_true")
    pa.add_argument("--normalize",
                    choices=["minmax", "vector", "area", "amide_I"],
                    default="amide_I")
    pa.add_argument("--deconv", action="store_true",
                    help="also run amide I deconvolution")
    pa.add_argument("--deconv-export", metavar="OUT.csv",
                    help="export full amide I deconvolution data to CSV "
                         "(wavenumber, observed, fit, 6 sub-bands, residual) "
                         "for plotting")
    pa.add_argument("--galalithe", action="store_true",
                    help="compute galalithe polymerization indicators")
    pa.add_argument("--reference", "-r", metavar="FILE",
                    help="path to a reference spectrum (theoretical or "
                         "commercial) for comparison")
    pa.add_argument("--use-default-galalithe-ref", action="store_true",
                    help="compare with the bundled exploratory galalithe DFT fragment "
                         "reference (not a certified material standard)")
    pa.add_argument("--no-preprocess-ref", action="store_true",
                    help="skip preprocessing on the reference (use if it "
                         "is already a clean theoretical/normalised spectrum)")
    pa.add_argument("--ref-threshold", type=float, default=0.10,
                    help="relative area diff threshold (fraction) for band "
                         "match/excess/deficit flag (default 0.10)")
    pa.add_argument("--subtract-ref", metavar="OUT.csv",
                    help="also save (sample - scale*reference) as a CSV")
    pa.add_argument("--subtract-scale", type=float, default=1.0,
                    help="scale factor applied to reference before subtraction")
    pa.add_argument("--report", help="output report file (.txt or .html)")
    pa.set_defaults(func=cmd_analyze)

    # compare
    pc = sub.add_parser("compare", help="before/after spectrum comparison")
    pc.add_argument("before", help="path to 'before' spectrum")
    pc.add_argument("after",  help="path to 'after'  spectrum")
    pc.add_argument("--no-preprocess", action="store_true")
    pc.add_argument("--threshold", type=float, default=0.05,
                    help="relative change threshold (default 0.05)")
    pc.add_argument("--diff-out", help="export difference spectrum to CSV")
    pc.set_defaults(func=cmd_compare)

    # deconv
    pd = sub.add_parser("deconv", help="amide I band deconvolution")
    pd.add_argument("input")
    pd.add_argument("--wn-min", type=float, default=1600.0)
    pd.add_argument("--wn-max", type=float, default=1700.0)
    pd.add_argument("--no-preprocess", action="store_true")
    pd.add_argument("--export", "-o", metavar="OUT.csv",
                    help="export deconvolution data (wavenumber, observed, "
                         "fit, 6 sub-bands, residual) to a CSV file for plotting")
    pd.set_defaults(func=cmd_deconv)

    # simulate
    ps = sub.add_parser("simulate", help="generate a synthetic spectrum")
    ps.add_argument("--kind", choices=["casein", "galalithe"],
                    default="casein")
    ps.add_argument("--species", choices=["bovine", "caprine", "ovine"],
                    default="bovine")
    ps.add_argument("--wn-min", type=float, default=600.0)
    ps.add_argument("--wn-max", type=float, default=4000.0)
    ps.add_argument("--crosslink", type=float, default=0.5,
                    help="synthetic marker amplitude parameter [0,1] (galalithe only)")
    ps.add_argument("--formol", type=float, default=0.05,
                    help="synthetic carbonyl amplitude parameter [0,1], not a measured concentration")
    ps.add_argument("--noise", type=float, default=0.0)
    ps.add_argument("--output", required=True)
    ps.set_defaults(func=cmd_simulate)

    # info
    pi = sub.add_parser("info", help="list known bands and references")
    pi.set_defaults(func=cmd_info)

    # gaussian
    pg = sub.add_parser("gaussian",
                          help="convert/inspect Gaussian .log frequency output")
    pg.add_argument("action",
                    choices=["convert", "inspect", "list-scales"],
                    help="'convert' = .log -> CSV spectrum, "
                         "'inspect' = print strongest modes, "
                         "'list-scales' = print known scale factors")
    pg.add_argument("input", nargs="?", help="Gaussian .log/.out file")
    pg.add_argument("--scale", type=float, default=None,
                    help="empirical scale factor (auto if omitted)")
    pg.add_argument("--level", default=None,
                    help="override detected DFT level (e.g. 'B3LYP/6-31G(d)')")
    pg.add_argument("--wn-min", type=float, default=400.0)
    pg.add_argument("--wn-max", type=float, default=4000.0)
    pg.add_argument("--resolution", type=float, default=2.0)
    pg.add_argument("--fwhm", type=float, default=15.0,
                    help="FWHM of each mode after convolution (cm-1)")
    pg.add_argument("--lineshape", choices=["gaussian", "lorentzian"],
                    default="lorentzian")
    pg.add_argument("--top-n", type=int, default=20,
                    help="for 'inspect' : number of strongest modes to print")
    pg.add_argument("--output", "-o",
                    help="for 'convert' : output CSV path")
    pg.set_defaults(func=cmd_gaussian)

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "analyze":
        if args.reference and args.use_default_galalithe_ref:
            parser.error("choose --reference OR --use-default-galalithe-ref")
        if args.deconv_export and not args.deconv:
            parser.error("--deconv-export requires --deconv")
        if args.subtract_ref and not (args.reference or args.use_default_galalithe_ref):
            parser.error("--subtract-ref requires a reference")
    if args.command == "gaussian":
        if args.action != "list-scales" and not args.input:
            parser.error("gaussian inspect/convert requires an input file")
        if args.action == "convert" and not args.output:
            parser.error("gaussian convert requires --output")
    try:
        return args.func(args)
    except (OSError, ValueError, ImportError, RuntimeError) as error:
        parser.exit(2, f"casein-ftir: error: {error}\n")


if __name__ == "__main__":
    main()
