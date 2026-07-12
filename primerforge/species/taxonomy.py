"""
NCBI taxonomy resolver.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio import Entrez


@dataclass(slots=True)
class TaxonomyRecord:

    taxid: str

    scientific_name: str

    rank: str

    lineage: list[str]

    lineage_ids: list[str]

    lineage_ranks: list[str]


class TaxonomyResolver:

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
    ) -> str:

        handle = Entrez.esearch(

            db="taxonomy",

            term=species,

            retmax=1,

        )

        record = Entrez.read(handle)

        handle.close()

        ids = record["IdList"]

        if not ids:

            raise ValueError(
                f"Species not found: {species}"
            )

        return ids[0]

    def fetch(
        self,
        taxid: str,
    ) -> TaxonomyRecord:

        handle = Entrez.efetch(

            db="taxonomy",

            id=taxid,

            retmode="xml",

        )

        records = Entrez.read(handle)

        handle.close()

        record = records[0]

        lineage = []
        lineage_ids = []
        lineage_ranks = []

        for item in record["LineageEx"]:

            lineage.append(
                item["ScientificName"]
            )

            lineage_ids.append(
                item["TaxId"]
            )

            lineage_ranks.append(
                item["Rank"]
            )

        return TaxonomyRecord(

            taxid=taxid,

            scientific_name=record["ScientificName"],

            rank=record["Rank"],

            lineage=lineage,

            lineage_ids=lineage_ids,

            lineage_ranks=lineage_ranks,

        )

    def resolve(
        self,
        species: str,
    ) -> TaxonomyRecord:

        taxid = self.search(
            species,
        )

        return self.fetch(
            taxid,
        )
