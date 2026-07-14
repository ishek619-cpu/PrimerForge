"""
Species workflow.

Coordinates the complete species-specific preprocessing
workflow before primer discovery.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from envoprimer.species.alignment import (
    AlignmentBuilder,
    AlignmentResult,
)
from envoprimer.species.collector import (
    SpeciesCollector,
    SpeciesDataset,
)
from envoprimer.species.diagnostics import (
    DiagnosticFinder,
    DiagnosticSite,
)
from envoprimer.species.windows import (
    DiagnosticWindow,
    DiagnosticWindowBuilder,
)


@dataclass(slots=True)
class WorkflowResult:
    """
    Complete preprocessing result.
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

        ####################################################
        # DATA COLLECTION
        ####################################################

        dataset = self.collector.collect(

            species=species,

            marker=marker,

            output_directory=output_directory,

        )

        ####################################################
        # ALIGNMENT
        ####################################################

        alignment = self.alignment.build(

            target_fasta=dataset.target_fasta,

            background_fasta=dataset.background_fasta,

            output_directory=output_directory,

            marker=marker,

        )

        ####################################################
        # DIAGNOSTIC SNP DISCOVERY
        ####################################################

        sites = self.diagnostics.find(

            alignment.target_alignment,

            alignment.background_alignment,

        )

        ####################################################
        # REFERENCE
        ####################################################

        record = next(

            SeqIO.parse(

                dataset.target_fasta,

                "fasta",

            )

        )

        ####################################################
        # WINDOWS
        ####################################################

        windows = self.windows.build(

            reference_sequence=str(

                record.seq,

            ),

            diagnostic_sites=sites,

        )

        ####################################################
        # RETURN
        ####################################################

        return WorkflowResult(

            dataset=dataset,

            alignment=alignment,

            diagnostic_sites=sites,

            diagnostic_windows=windows,

        )
