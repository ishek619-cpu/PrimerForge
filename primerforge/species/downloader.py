"""
NCBI species sequence downloader.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio import Entrez
from Bio import SeqIO

from primerforge.species.cache import DownloadCache
from primerforge.species.markers import build_query


@dataclass(slots=True)
class DownloadResult:

    species: str

    marker: str

    accessions: list[str]

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

        return list(record["IdList"])

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

        with output.open(
            "w",
            encoding="utf-8",
        ) as out_handle:

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

                SeqIO.write(
                    SeqIO.parse(
                        handle,
                        "fasta",
                    ),
                    out_handle,
                    "fasta",
                )

                handle.close()

        return output

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

        #
        # Use cached FASTA if available.
        #
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

        return DownloadResult(

            species=species,

            marker=marker,

            accessions=ids,

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

                results.append(result)

            except Exception as exc:

                print(
                    f"Skipping {species}: {exc}"
                )

        return results

    def fasta_files(
        self,
        downloads: list[DownloadResult],
    ) -> list[Path]:
        """
        Extract FASTA paths from download results.
        """

        return [
            result.output_fasta
            for result in downloads
        ]
