"""
Species collection workflow.

Phase 1:
- Resolve taxonomy
- Find related species
- Download target species
- Download related species
- Merge background FASTA
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from primerforge.io.merge import FASTAMerger
from primerforge.species.downloader import (
    DownloadResult,
    SpeciesDownloader,
)
from primerforge.species.relatives import (
    RelativeSpeciesFinder,
)
from primerforge.species.taxonomy import (
    TaxonomyResolver,
)


@dataclass(slots=True)
class SpeciesDataset:
    """
    Species dataset produced by the collector.
    """

    species: str

    marker: str

    taxonomy: object

    genus: str

    target_result: DownloadResult

    background_results: list[DownloadResult]

    target_fasta: Path

    background_fasta: Path


class SpeciesCollector:
    """
    Collect all sequence data required for primer design.
    """

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
    ):

        self.taxonomy = TaxonomyResolver(
            email=email,
            api_key=api_key,
        )

        self.relatives = RelativeSpeciesFinder(
            email=email,
            api_key=api_key,
        )

        self.downloader = SpeciesDownloader(
            email=email,
            api_key=api_key,
        )

        self.merger = FASTAMerger()

    def collect(
        self,
        species: str,
        marker: str,
        output_directory: Path,
    ) -> SpeciesDataset:

        output_directory.mkdir(
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
        # Determine genus
        #
        genus = self.relatives.genus(
            taxonomy,
        )

        #
        # Find related species
        #
        relatives = self.relatives.find(
            taxonomy,
        )

        print(
            f"Found {len(relatives)} related species."
        )

        #
        # Download target
        #
        print(
            f"Downloading target: {species}"
        )

        target = self.downloader.download(
            species=species,
            marker=marker,
            output_directory=output_directory,
        )

        #
        # Download background species
        #
        print(
            "Downloading related species..."
        )

        background = self.downloader.download_many(
            species_list=relatives,
            marker=marker,
            output_directory=output_directory,
        )

        #
        # Merge all background FASTA files
        #
        background_fasta = (
            output_directory
            / "background.fasta"
        )

        self.merger.merge(
            self.downloader.fasta_files(
                background,
            ),
            background_fasta,
        )

        return SpeciesDataset(

            species=species,

            marker=marker,

            taxonomy=taxonomy,

            genus=genus,

            target_result=target,

            background_results=background,

            target_fasta=target.output_fasta,

            background_fasta=background_fasta,

        )
