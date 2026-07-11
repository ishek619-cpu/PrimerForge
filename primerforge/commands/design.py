"""
PrimerForge design command.
"""

from __future__ import annotations

from pathlib import Path

from primerforge.primer.discovery import PrimerDiscovery
from primerforge.reference.reference import Reference

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

    print(f"Reference : {reference_fasta}")

    if alignment is not None:
        print(f"Alignment : {alignment}")

    if annotation is not None:
        print(f"Annotation: {annotation}")

    if gene:
        print(f"Gene       : {gene}")

    print(f"Output     : {output}")

    print()

    #
    # Build Reference object
    #
    reference = Reference(
        fasta=reference_fasta,
        alignment=alignment,
        annotation=annotation,
    )

    #
    # Show annotation summary
    #
    if reference.has_annotation:

        print(
            f"Loaded {len(reference.genes)} annotated features."
        )

        if gene:

            feature = reference.get_gene(
                gene,
            )

            if feature is None:

                print()

                print(
                    f"ERROR: Gene '{gene}' was not found."
                )

                print()

                print(
                    "Available genes:"
                )

                for name in reference.gene_names():

                    print(
                        f"  - {name}"
                    )

                return 1

            print(
                f"Target gene : {feature.name}"
            )

            print(
                f"Coordinates : {feature.start}-{feature.end}"
            )

    print()

    print("Designing primers...")

    discovery = PrimerDiscovery()

    #
    # NOTE:
    # Gene-aware discovery will be implemented next.
    # For now the Reference object is already available.
    #
    pairs = discovery.discover(
        reference_fasta=reference.fasta,
        alignment_fasta=reference.alignment,
    )

    print(
        f"Found {len(pairs)} primer pairs."
    )

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

    print(
        f"  HTML  : {output/'primerforge_report.html'}"
    )

    print(
        f"  CSV   : {output/'primerforge_report.csv'}"
    )

    print(
        f"  JSON  : {output/'primerforge_report.json'}"
    )

    print(
        f"  Excel : {output/'primerforge_report.xlsx'}"
    )

    print()

    return 0
