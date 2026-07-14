"""
NCBI species sequence downloader.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio import Entrez
from Bio import SeqIO

from envoprimer.species.cache import DownloadCache
from envoprimer.species.filter import SequenceFilter
from envoprimer.species.markers import build_query
from envoprimer.species.metadata import SequenceMetadata


@dataclass(slots=True)
class DownloadResult:

    species: str

    marker: str

    accessions: list[str]

    metadata: list[SequenceMetadata]

    output_fasta: Path


class SpeciesDownloader:

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
    ):

        Entrez.email = email

        if api_key:
            Entrez.api_key = api_key

        self.cache = DownloadCache(
            Path("cache"),
        )

        self.filter = SequenceFilter(
            minimum_length=500,
        )

    def search(
        self,
        species: str,
        marker: str,
    ) -> list[str]:

        query = build_query(
            species,
            marker,
        )

        handle = Entrez.esearch(
            db="nucleotide",
            term=query,
            retmax=100000,
        )

        record = Entrez.read(handle)

        handle.close()

        ids = list(record["IdList"])

        print(
            f"Found {len(ids)} NCBI records."
        )

        return ids

    def fetch(
        self,
        accessions: Iterable[str],
        output: Path,
        batch_size: int = 500,
    ) -> Path:

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        ids = list(accessions)

        if not ids:

            output.write_text(
                "",
                encoding="utf-8",
            )

            return output

        records = []

        for start in range(
            0,
            len(ids),
            batch_size,
        ):

            batch = ids[
                start:start + batch_size
            ]

            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(batch),
                rettype="fasta",
                retmode="text",
            )

            records.extend(

                list(

                    SeqIO.parse(
                        handle,
                        "fasta",
                    )

                )

            )

            handle.close()

        print(
            f"Downloaded {len(records)} sequences."
        )

        records = self.filter.clean(
            records,
        )

        print(
            f"Retained {len(records)} high-quality sequences."
        )

        with output.open(
            "w",
            encoding="utf-8",
        ) as handle:

            SeqIO.write(
                records,
                handle,
                "fasta",
            )

        return output

    def build_metadata(
        self,
        fasta: Path,
        species: str,
        marker: str,
        role: str,
    ) -> list[SequenceMetadata]:

        metadata = []

        for record in SeqIO.parse(
            fasta,
            "fasta",
        ):

            accession = record.id.split(".")[0]

            metadata.append(

                SequenceMetadata(

                    accession=accession,

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
    ) -> DownloadResult:

        cache_file = self.cache.path(
            species,
            marker,
        )

        if cache_file.exists():

            print(
                f"Using cached sequences for {species}"
            )

            ids = []

        else:

            ids = self.search(
                species,
                marker,
            )

            self.fetch(
                ids,
                cache_file,
            )

        directory = (
            output_directory
            / species.replace(
                " ",
                "_",
            )
        )

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        fasta = directory / "target.fasta"

        fasta.write_bytes(
            cache_file.read_bytes()
        )

        metadata = self.build_metadata(
            fasta=fasta,
            species=species,
            marker=marker,
            role="target",
        )

        print(
            f"Final target dataset: {len(metadata)} sequences."
        )

        return DownloadResult(

            species=species,

            marker=marker,

            accessions=ids,

            metadata=metadata,

            output_fasta=fasta,

        )

    def download_many(
        self,
        species_list: list[str],
        marker: str,
        output_directory: Path,
    ) -> list[DownloadResult]:

        results = []

        for species in species_list:

            try:

                result = self.download(
                    species=species,
                    marker=marker,
                    output_directory=output_directory,
                )

                for record in result.metadata:

                    record.role = "background"

                results.append(
                    result,
                )

            except Exception as exc:

                print(
                    f"Skipping {species}: {exc}"
                )

        return results

    def fasta_files(
        self,
        downloads: list[DownloadResult],
    ) -> list[Path]:

        return [

            result.output_fasta

            for result in downloads

        ]
