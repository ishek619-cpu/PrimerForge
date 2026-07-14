"""
EnvoPrimer design command.
"""

from __future__ import annotations

from pathlib import Path

from envoprimer.core.pipeline import Pipeline
from envoprimer.core.species_pipeline import SpeciesPipeline

from envoprimer.primer.discovery import PrimerDiscovery
from envoprimer.reference.reference import Reference

from envoprimer.report.csv import CSVReport
from envoprimer.report.excel import ExcelReport
from envoprimer.report.html import HTMLReport
from envoprimer.report.json import JSONReport


def run_design(args) -> int:

    #
    # Species workflow
    #
    if getattr(args, "config", None):

        pipeline = SpeciesPipeline()

        pipeline.run(
            Path(args.config),
        )

        return 0

    #
    # Legacy FASTA workflow
    #

    print()

    print("EnvoPrimer")
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
