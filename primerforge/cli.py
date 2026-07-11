"""
PrimerForge command-line interface.
"""

from __future__ import annotations

import argparse
import sys

from primerforge.commands.design import run_design
from primerforge.commands.report import run_report
from primerforge.commands.validate import run_validate


def build_parser() -> argparse.ArgumentParser:
    """
    Build the command-line parser.
    """

    parser = argparse.ArgumentParser(
        prog="primerforge",
        description="PrimerForge: Advanced PCR primer design and validation platform.",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    ####################################################################
    # DESIGN
    ####################################################################

    design = subparsers.add_parser(
        "design",
        help="Design PCR primer pairs.",
    )

    design.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA file.",
    )

    design.add_argument(
        "--alignment",
        help="Population alignment FASTA.",
    )

    design.add_argument(
        "--annotation",
        help="GenBank annotation (.gb/.gbk).",
    )

    design.add_argument(
        "--gene",
        help="Design primers for a specific annotated gene (e.g. COI, CYTB, 12S).",
    )

    design.add_argument(
        "--output",
        default="results",
        help="Output directory.",
    )

    design.set_defaults(
        func=run_design,
    )

    ####################################################################
    # VALIDATE
    ####################################################################

    validate = subparsers.add_parser(
        "validate",
        help="Validate primer pairs.",
    )

    validate.add_argument(
        "--input",
        required=True,
        help="Primer input file.",
    )

    validate.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA.",
    )

    validate.add_argument(
        "--alignment",
        help="Population alignment.",
    )

    validate.add_argument(
        "--annotation",
        help="GenBank annotation.",
    )

    validate.add_argument(
        "--output",
        default="results",
        help="Output directory.",
    )

    validate.set_defaults(
        func=run_validate,
    )

    ####################################################################
    # REPORT
    ####################################################################

    report = subparsers.add_parser(
        "report",
        help="Generate reports.",
    )

    report.add_argument(
        "--input",
        required=True,
        help="Input results.",
    )

    report.add_argument(
        "--output",
        default="results",
        help="Output directory.",
    )

    report.set_defaults(
        func=run_report,
    )

    return parser


def main() -> int:
    """
    CLI entry point.
    """

    parser = build_parser()

    args = parser.parse_args()

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
