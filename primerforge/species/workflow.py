"""
Species workflow.

Coordinates the complete species-specific preprocessing
workflow before primer discovery.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from primerforge.species.alignment import (
    AlignmentBuilder,
    AlignmentResult,
)
from primerforge.species.collector import (
    SpeciesCollector,
    SpeciesDataset,
)
from primerforge.species.diagnostics import (
    DiagnosticFinder,
    DiagnosticSite,
)
from primerforge.species.windows import (
    DiagnosticWindow,
    DiagnosticWindowBuilder,
)


@dataclass(slots=True)
class WorkflowResult:
    """
    Complete preprocessing result for primer discovery.
    """

    dataset: SpeciesDataset

    alignment: AlignmentResult

    diagnostic_sites: list[DiagnosticSite]

    diagnostic_windows: list[DiagnosticWindow]


class SpeciesWorkflow:
    """
    Complete preprocessing workflow.
    """

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
    ):

        self.collector = SpeciesCollector(
            email=email,
            api_key=api_key,
        )

        self.alignment = AlignmentBuilder()

        self.diagnostics = DiagnosticFinder()

        self.windows = DiagnosticWindowBuilder()

    def run(
        self,
        species: str,
        marker: str,
        output_directory: Path,
    ) -> WorkflowResult:

        #
        # Collect sequences
        #
        dataset = self.collector.collect(
            species=species,
            marker=marker,
            output_directory=output_directory,
        )

        #
        # Build alignments
        #
        alignment = self.alignment.build(
            target_fasta=dataset.target_fasta,
            background_fasta=dataset.background_fasta,
            output_directory=output_directory,
        )

        #
        # Diagnostic SNPs
        #
        sites = self.diagnostics.find(
            alignment.target_alignment,
            alignment.background_alignment,
        )

        #
        # Load reference sequence
        #
        record = next(
            SeqIO.parse(
                dataset.target_fasta,
                "fasta",
            )
        )

        #
        # Diagnostic windows
        #
        windows = self.windows.build(
            reference_sequence=str(record.seq),
            diagnostic_sites=sites,
        )

        return WorkflowResult(

            dataset=dataset,

            alignment=alignment,

            diagnostic_sites=sites,

            diagnostic_windows=windows,

        )
