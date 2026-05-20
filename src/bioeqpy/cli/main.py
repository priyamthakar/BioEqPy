"""Command-line interface for BioEqPy."""

from __future__ import annotations

import argparse

from bioeqpy import __version__
from bioeqpy.analysis import analyze
from bioeqpy.power import achieved_power, sample_size
from bioeqpy.reporting import format_ci_table


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level CLI parser."""
    parser = argparse.ArgumentParser(prog="bioeqpy")
    parser.add_argument("--version", action="store_true", help="Show BioEqPy version.")
    subparsers = parser.add_subparsers(dest="command")

    analyze_parser = subparsers.add_parser("analyze", help="Run BE analysis on a CSV or Excel file.")
    analyze_parser.add_argument("--data", required=True, help="Input CSV or Excel file.")
    analyze_parser.add_argument("--parameters", nargs="*", default=None, help="PK parameter columns to analyze.")
    analyze_parser.add_argument("--report", default=None, help="Optional Markdown report output path.")

    power_parser = subparsers.add_parser("power", help="Estimate 2x2 crossover sample size.")
    power_parser.add_argument("--design", default="2x2")
    power_parser.add_argument("--cv", type=float, required=True)
    power_parser.add_argument("--gmr", type=float, default=0.95)
    power_parser.add_argument("--power", type=float, default=0.80)
    power_parser.add_argument("--dropout-rate", type=float, default=0.0)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the BioEqPy CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version or args.command == "version":
        print(f"bioeqpy {__version__}")
        return 0

    if args.command == "analyze":
        results = analyze(args.data, parameters=args.parameters, report=args.report)
        print(format_ci_table(results).to_string(index=False))
        return 0

    if args.command == "power":
        n = sample_size(
            design=args.design,
            cv=args.cv,
            gmr=args.gmr,
            power=args.power,
            dropout_rate=args.dropout_rate,
        )
        observed_power = achieved_power(n=n, cv=args.cv, gmr=args.gmr)
        print(f"Required N: {n}")
        print(f"Approximate power at N={n}: {observed_power:.4f}")
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

