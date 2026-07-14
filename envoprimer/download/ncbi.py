"""
NCBI client for EnvoPrimer.

All communication with Entrez happens here.
"""

from __future__ import annotations

from typing import Iterable

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class NCBIClient:

    def __init__(

        self,

        email: str,

        api_key: str | None = None,

    ):

        Entrez.email = email

        if api_key:

            Entrez.api_key = api_key

    ####################################################################
    # Search
    ####################################################################

    def search(

        self,

        query: str,

        database: str = "nucleotide",

        maximum: int = 100000,

    ) -> list[str]:

        handle = Entrez.esearch(

            db=database,

            term=query,

            retmax=maximum,

        )

        result = Entrez.read(handle)

        handle.close()

        return list(result["IdList"])

    ####################################################################
    # FASTA download
    ####################################################################

    def fetch_fasta(

        self,

        ids: Iterable[str],

        batch_size: int = 500,

    ) -> list[SeqRecord]:

        ids = list(ids)

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

        return records

    ####################################################################
    # GenBank download
    ####################################################################

    def fetch_genbank(

        self,

        ids: Iterable[str],

        batch_size: int = 200,

    ) -> list[SeqRecord]:

        ids = list(ids)

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

                rettype="gb",

                retmode="text",

            )

            records.extend(

                list(

                    SeqIO.parse(

                        handle,

                        "genbank",

                    )

                )

            )

            handle.close()

        return records
