"""
EnvoPrimer command-line interface.
"""

from __future__ import annotations

import argparse
import sys

from envoprimer.commands.design import run_design
from envoprimer.commands.report import run_report
from envoprimer.commands.validate import run_validate


def build_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(
        prog="envoprimer",
        description="EnvoPrimer: Species-specific eDNA primer design platform.",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    ############################################################
    # DESIGN
    ############################################################

    design = subparsers.add_parser(
        "design",
        help="Design species-specific eDNA primers.",
    )

    #
    # Reference workflow
    #
    design.add_argument(
        "--reference",
        help="Reference FASTA file.",
    )

    design.add_argument(
        "--alignment",
        help="Population alignment FASTA.",
    )

    design.add_argument(
        "--annotation",
        help="Reference GenBank annotation.",
    )

    design.add_argument(
        "--gene",
        help="Target gene.",
    )

    #
    # Species workflow
    #
    design.add_argument(
        "--config",
        help="EnvoPrimer YAML configuration.",
    )

    design.add_argument(
        "--output",
        default="results",
        help="Output directory.",
    )

    design.set_defaults(
        func=run_design,
    )

    ############################################################
    # VALIDATE
    ############################################################

    validate = subparsers.add_parser(
        "validate",
        help="Validate designed primers.",
    )

    validate.add_argument(
        "--input",
        required=True,
        help="Primer file.",
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

    ############################################################
    # REPORT
    ############################################################

    report = subparsers.add_parser(
        "report",
        help="Generate EnvoPrimer reports.",
    )

    report.add_argument(
        "--input",
        required=True,
        help="Input JSON file.",
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


def main():

    parser = build_parser()

    args = parser.parse_args()

    return args.func(args)


if __name__ == "__main__":

    sys.exit(main())
