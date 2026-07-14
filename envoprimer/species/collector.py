"""
Species collection workflow.

Collects, cleans and prepares all datasets required for
species-specific primer design.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from envoprimer.download.models import DownloadResult
from envoprimer.download.pipeline import DownloadPipeline
from envoprimer.io.merge import FASTAMerger
from envoprimer.species.cleaner import SequenceCleaner
from envoprimer.species.metadata import SequenceMetadata
from envoprimer.species.relatives import RelativeSpeciesFinder
from envoprimer.species.taxonomy import TaxonomyResolver


@dataclass(slots=True)
class SpeciesDataset:
    """
    Species dataset produced by the collector.
    """

    species: str

    marker: str

    taxonomy: object

    genus: tuple[str, str]

    target_result: DownloadResult

    background_results: list[DownloadResult]

    target_metadata: list[SequenceMetadata]

    background_metadata: list[SequenceMetadata]

    target_fasta: Path

    background_fasta: Path


class SpeciesCollector:
    """
    Collect all sequence data required for primer discovery.
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

        self.downloader = DownloadPipeline(
            email=email,
            api_key=api_key,
        )

        self.cleaner = SequenceCleaner()

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

        taxonomy = self.taxonomy.resolve(
            species,
        )

        genus = self.relatives.genus(
            taxonomy,
        )

        relatives = self.relatives.find(
            taxonomy,
        )

        print(
            f"Found {len(relatives)} related species."
        )

        print(
            f"Downloading target: {species}"
        )

        target = self.downloader.download(
            species=species,
            marker=marker,
            output_directory=output_directory,
            role="target",
        )

        target_clean = (
            output_directory
            / species.replace(" ", "_")
            / "target.cleaned.fasta"
        )

        kept = self.cleaner.clean(
            target.output_fasta,
            target_clean,
        )

        print(
            f"Target sequences kept: {kept}"
        )

        target.output_fasta = target_clean

        print(
            "Downloading related species..."
        )

        background = self.downloader.download_many(
            species_list=relatives,
            marker=marker,
            output_directory=output_directory,
        )

        cleaned_fastas = []

        background_metadata = []

        for result in background:

            cleaned = result.output_fasta.with_suffix(
                ".cleaned.fasta"
            )

            kept = self.cleaner.clean(
                result.output_fasta,
                cleaned,
            )

            if kept == 0:
                continue

            result.output_fasta = cleaned

            cleaned_fastas.append(
                cleaned
            )

            background_metadata.extend(
                result.metadata
            )

        background_fasta = (
            output_directory
            / "background.cleaned.fasta"
        )

        self.merger.merge(
            cleaned_fastas,
            background_fasta,
        )

        return SpeciesDataset(

            species=species,

            marker=marker,

            taxonomy=taxonomy,

            genus=genus,

            target_result=target,

            background_results=background,

            target_metadata=target.metadata,

            background_metadata=background_metadata,

            target_fasta=target.output_fasta,

            background_fasta=background_fasta,

        )
