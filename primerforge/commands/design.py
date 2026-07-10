"""
PrimerForge design command.
"""

from __future__ import annotations

from pathlib import Path

from primerforge.primer.discovery import PrimerDiscovery

from primerforge.report.csv import CSVReport
from primerforge.report.excel import ExcelReport
from primerforge.report.html import HTMLReport
from primerforge.report.json import JSONReport


def run_design(args) -> int:
    """
    Run the complete primer discovery workflow.
    """

    print()

    print("PrimerForge")
    print("=" * 60)

    reference = Path(args.reference)

    alignment = (
        Path(args.alignment)
        if args.alignment
        else None
    )

    output = Path(args.output)

    output.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(f"Reference : {reference}")

    if alignment is not None:
        print(f"Alignment : {alignment}")

    print(f"Output     : {output}")

    print()

    print("Designing primers...")

    discovery = PrimerDiscovery()

    pairs = discovery.discover(
        reference_fasta=reference,
        alignment_fasta=alignment,
    )

    print(f"Found {len(pairs)} primer pairs.")

    print()

    print("Writing reports...")

    HTMLReport().write(
        pairs,
        output / "primerforge_report.html",
    )

    CSVReport().write(
        pairs,
        output / "primerforge_report.csv",
    )

    JSONReport().write(
        pairs,
        output / "primerforge_report.json",
    )

    ExcelReport().write(
        pairs,
        output / "primerforge_report.xlsx",
    )

    print("Done.")

    if pairs:

        print()

        print(
            f"Best score : {pairs[0].score:.2f}"
        )

        print(
            f"Best pair  : {pairs[0].forward.sequence}"
        )

        print(
            f"             {pairs[0].reverse.sequence}"
        )

    print()

    print("Reports")

    print(f"  HTML  : {output/'primerforge_report.html'}")
    print(f"  CSV   : {output/'primerforge_report.csv'}")
    print(f"  JSON  : {output/'primerforge_report.json'}")
    print(f"  Excel : {output/'primerforge_report.xlsx'}")

    print()

    return 0
