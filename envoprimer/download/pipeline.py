"""
EnvoPrimer download pipeline.

Downloads, extracts and prepares marker datasets.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from envoprimer.download.downloader import SequenceDownloader
from envoprimer.download.models import DownloadResult
from envoprimer.extraction.gene import GeneExtractor
from envoprimer.species.metadata import SequenceMetadata


class DownloadPipeline:

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
    ):

        self.downloader = SequenceDownloader(
            email=email,
            api_key=api_key,
        )

        self.extractor = GeneExtractor()

    def _metadata(
        self,
        species: str,
        marker: str,
        fasta: Path,
        role: str,
    ) -> list[SequenceMetadata]:

        metadata = []

        for record in SeqIO.parse(
            fasta,
            "fasta",
        ):

            metadata.append(

                SequenceMetadata(

                    accession=record.id.split(".")[0],

                    species=species,

                    marker=marker,

                    definition=record.description,

                    length=len(record.seq),

                    source="NCBI",

                    role=role,

                    fasta=str(fasta),

                )

            )

        return metadata

    def download(
        self,
        species: str,
        marker: str,
        output_directory: Path,
        role: str = "target",
    ) -> DownloadResult:

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        fasta = (
            output_directory
            / species.replace(" ", "_")
            / "target.fasta"
        )

        records = self.downloader.download_genbank(
            species,
            marker,
        )

        gb = fasta.with_suffix(".gb")

        gb.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        SeqIO.write(
            records,
            gb,
            "genbank",
        )

        genes = self.extractor.extract(
            gb,
            marker,
        )

        self.extractor.write_fasta(
            genes,
            fasta,
        )

        metadata = self._metadata(
            species,
            marker,
            fasta,
            role,
        )

        return DownloadResult(

            species=species,

            marker=marker,

            accessions=[
                g.accession
                for g in genes
            ],

            metadata=metadata,

            output_fasta=fasta,

        )

    def download_many(
        self,
        species_list: list[str],
        marker: str,
        output_directory: Path,
    ) -> list[DownloadResult]:
        """
        Download multiple background species.
        """

        results: list[DownloadResult] = []

        print()
        print("=" * 60)
        print("DOWNLOADING BACKGROUND SPECIES")
        print("=" * 60)

        for species in species_list:

            try:

                print(f"Downloading: {species}")

                result = self.download(
                    species=species,
                    marker=marker,
                    output_directory=output_directory,
                    role="background",
                )

                results.append(result)

            except Exception as exc:

                print(
                    f"Skipping {species}: {exc}"
                )

        print()
        print(
            f"Successfully downloaded {len(results)} background species."
        )

        return results
