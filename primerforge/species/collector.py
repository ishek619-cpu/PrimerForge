"""
Species collection workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from primerforge.species.aligner import MAFFTAligner
from primerforge.species.cleaner import SequenceCleaner
from primerforge.species.diagnostics import DiagnosticFinder
from primerforge.species.downloader import SpeciesDownloader
from primerforge.species.markers import MarkerDatabase
from primerforge.species.relatives import RelativeFinder
from primerforge.species.taxonomy import TaxonomyResolver


@dataclass(slots=True)
class SpeciesDataset:
    """
    Complete dataset used for species-specific primer design.
    """

    species: str

    marker: str

    taxonomy: object

    relatives: list

    target_fasta: Path

    relative_fastas: list[Path]

    cleaned_target: Path

    cleaned_relatives: list[Path]

    alignment: Path

    diagnostics: list


class SpeciesCollector:
    """
    Complete workflow for collecting all data required for
    species-specific primer design.
    """

    def __init__(self):

        self.taxonomy = TaxonomyResolver()

        self.relatives = RelativeFinder()

        self.markers = MarkerDatabase()

        self.downloader = SpeciesDownloader()

        self.cleaner = SequenceCleaner()

        self.aligner = MAFFTAligner()

        self.diagnostics = DiagnosticFinder()

    def collect(
        self,
        species: str,
        marker: str,
        output: Path,
    ) -> SpeciesDataset:

        output.mkdir(
            parents=True,
            exist_ok=True,
        )

        #
        # Resolve taxonomy
        #
        taxonomy = self.taxonomy.resolve(
            species,
        )

        #
        # Find relatives
        #
        relatives = self.relatives.find(
            taxonomy,
        )

        #
        # Download target sequences
        #
        target_fasta = self.downloader.download(
            species=species,
            marker=marker,
            output=output / "target.fasta",
        )

        #
        # Download relatives
        #
        relative_fastas = []

        for relative in relatives:

            fasta = self.downloader.download(
                species=relative,
                marker=marker,
                output=output / f"{relative.replace(' ','_')}.fasta",
            )

            relative_fastas.append(
                fasta,
            )

        #
        # Clean target
        #
        cleaned_target = self.cleaner.clean(
            target_fasta,
            output / "target.cleaned.fasta",
        )

        #
        # Clean relatives
        #
        cleaned_relatives = []

        for fasta in relative_fastas:

            cleaned = self.cleaner.clean(
                fasta,
                fasta.with_suffix(".cleaned.fasta"),
            )

            cleaned_relatives.append(
                cleaned,
            )

        #
        # Merge all cleaned FASTA files
        #
        merged = output / "merged.fasta"

        with merged.open(
            "w",
            encoding="utf-8",
        ) as out:

            out.write(
                cleaned_target.read_text(
                    encoding="utf-8",
                )
            )

            for fasta in cleaned_relatives:

                out.write(
                    fasta.read_text(
                        encoding="utf-8",
                    )
                )

        #
        # Alignment
        #
        alignment = output / "alignment.fasta"

        self.aligner.align(
            merged,
            alignment,
        )

        #
        # Diagnostic SNPs
        #
        diagnostics = self.diagnostics.find(
            alignment,
        )

        return SpeciesDataset(

            species=species,

            marker=marker,

            taxonomy=taxonomy,

            relatives=relatives,

            target_fasta=target_fasta,

            relative_fastas=relative_fastas,

            cleaned_target=cleaned_target,

            cleaned_relatives=cleaned_relatives,

            alignment=alignment,

            diagnostics=diagnostics,

        )
