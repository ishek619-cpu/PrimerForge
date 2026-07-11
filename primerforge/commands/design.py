"""
PrimerForge design command.
"""

from __future__ import annotations

from pathlib import Path

from primerforge.core.pipeline import Pipeline

from primerforge.primer.discovery import PrimerDiscovery
from primerforge.reference.reference import Reference

from primerforge.report.csv import CSVReport
from primerforge.report.excel import ExcelReport
from primerforge.report.html import HTMLReport
from primerforge.report.json import JSONReport


def run_design(args) -> int:

    #
    # NEW SPECIES WORKFLOW
    #
    if getattr(args, "config", None):

        pipeline = Pipeline()

        pipeline.run(
            Path(args.config),
        )

        return 0

    #
    # LEGACY FASTA WORKFLOW
    #

    print()

    print("PrimerForge")
    print("=" * 60)

    reference_fasta = Path(args.reference)

    alignment = (
        Path(args.alignment)
        if args.alignment
        else None
    )

    annotation = (
        Path(args.annotation)
        if getattr(args, "annotation", None)
        else None
    )

    gene = getattr(
        args,
        "gene",
        None,
    )

    output = Path(args.output)

    output.mkdir(
        parents=True,
        exist_ok=True,
    )

    reference = Reference(
        fasta=reference_fasta,
        alignment=alignment,
        annotation=annotation,
    )

    discovery = PrimerDiscovery()

    pairs = discovery.discover(
        reference_fasta=reference.fasta,
        alignment_fasta=reference.alignment,
    )

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

    return 0
