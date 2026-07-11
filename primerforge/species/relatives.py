"""
Find related species from NCBI taxonomy.
"""

from __future__ import annotations

from Bio import Entrez


class RelativeSpeciesFinder:

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
    ):

        Entrez.email = email

        if api_key:

            Entrez.api_key = api_key

    def genus(
        self,
        taxonomy_record,
    ) -> str:

        for item in reversed(
            taxonomy_record.lineage,
        ):

            if item != taxonomy_record.scientific_name:

                return item

        raise ValueError(
            "Unable to determine genus."
        )

    def search(
        self,
        genus: str,
    ) -> list[str]:

        handle = Entrez.esearch(

            db="taxonomy",

            term=f"{genus}[Genus]",

            retmax=1000,

        )

        record = Entrez.read(
            handle,
        )

        handle.close()

        return list(
            record["IdList"],
        )
