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

    parser = argparse.ArgumentParser(
        prog="primerforge",
        description="PrimerForge: Advanced PCR primer design platform.",
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
        help="Design PCR primers.",
    )

    #
    # Old workflow
    #
    design.add_argument(
        "--reference",
        help="Reference FASTA.",
    )

    design.add_argument(
        "--alignment",
        help="Population alignment.",
    )

    design.add_argument(
        "--annotation",
        help="GenBank annotation.",
    )

    design.add_argument(
        "--gene",
        help="Target annotated gene.",
    )

    #
    # New workflow
    #
    design.add_argument(
        "--config",
        help="PrimerForge YAML configuration.",
    )

    design.add_argument(
        "--output",
        default="results",
    )

    design.set_defaults(
        func=run_design,
    )

    ############################################################
    # VALIDATE
    ############################################################

    validate = subparsers.add_parser(
        "validate",
        help="Validate primers.",
    )

    validate.add_argument(
        "--input",
        required=True,
    )

    validate.add_argument(
        "--reference",
        required=True,
    )

    validate.add_argument(
        "--alignment",
    )

    validate.add_argument(
        "--annotation",
    )

    validate.add_argument(
        "--output",
        default="results",
    )

    validate.set_defaults(
        func=run_validate,
    )

    ############################################################
    # REPORT
    ############################################################

    report = subparsers.add_parser(
        "report",
    )

    report.add_argument(
        "--input",
        required=True,
    )

    report.add_argument(
        "--output",
        default="results",
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
