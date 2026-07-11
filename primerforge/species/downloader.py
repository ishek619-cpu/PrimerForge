"""
NCBI species sequence downloader.

Downloads all sequences for a target species and marker
using the NCBI Entrez E-utilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio import Entrez
from Bio import SeqIO

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

    def search(
        self,
        species: str,
        marker: str,
    ) -> list[str]:
        """
        Search GenBank for all nucleotide accessions
        matching a species and marker.
        """

        query = build_query(
            species,
            marker,
        )

        handle = Entrez.esearch(

            db="nucleotide",

            term=query,

            retmax=100000,

        )

        record = Entrez.read(
            handle,
        )

        handle.close()

        return list(
            record["IdList"],
        )

    def fetch(
        self,
        accessions: Iterable[str],
        output: Path,
        batch_size: int = 500,
    ) -> Path:
        """
        Download FASTA sequences in batches.
        """

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        ids = list(
            accessions,
        )

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

            for i in range(
                0,
                len(ids),
                batch_size,
            ):

                batch = ids[
                    i:i + batch_size
                ]

                handle = Entrez.efetch(

                    db="nucleotide",

                    id=",".join(batch),

                    rettype="fasta",

                    retmode="text",

                )

                records = list(
                    SeqIO.parse(
                        handle,
                        "fasta",
                    )
                )

                handle.close()

                SeqIO.write(
                    records,
                    out_handle,
                    "fasta",
                )

        return output

    def download(
        self,
        species: str,
        marker: str,
        output_directory: Path,
    ) -> DownloadResult:

        ids = self.search(
            species,
            marker,
        )

        directory = (
            output_directory
            /
            species.replace(
                " ",
                "_",
            )
        )

        fasta = (
            directory
            /
            "target.fasta"
        )

        self.fetch(
            ids,
            fasta,
        )

        return DownloadResult(

            species=species,

            marker=marker,

            accessions=ids,

            output_fasta=fasta,

        )
