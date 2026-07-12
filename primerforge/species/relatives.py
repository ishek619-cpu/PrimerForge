"""
Find related species from NCBI taxonomy.
"""

from __future__ import annotations

from Bio import Entrez


class RelativeSpeciesFinder:
    """
    Find closely related species from the NCBI taxonomy.
    """

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
    ) -> tuple[str, str]:

        for name, taxid, rank in zip(

            taxonomy_record.lineage,

            taxonomy_record.lineage_ids,

            taxonomy_record.lineage_ranks,

        ):

            if rank == "genus":

                return (
                    name,
                    taxid,
                )

        raise ValueError(
            "Genus not found."
        )

    def search(
        self,
        genus_taxid: str,
        exclude: str | None = None,
    ) -> list[str]:

        handle = Entrez.esearch(

            db="taxonomy",

            term=f"txid{genus_taxid}[Subtree]",

            retmax=5000,

        )

        search = Entrez.read(handle)

        handle.close()

        ids = search["IdList"]

        if not ids:

            return []

        handle = Entrez.efetch(

            db="taxonomy",

            id=",".join(ids),

            retmode="xml",

        )

        records = Entrez.read(handle)

        handle.close()

        species = []

        for record in records:

            if record["Rank"] != "species":

                continue

            name = record["ScientificName"]

            if exclude is not None and name == exclude:

                continue

            species.append(name)

        return sorted(set(species))

    def find(
        self,
        taxonomy_record,
    ) -> list[str]:

        genus_name, genus_taxid = self.genus(
            taxonomy_record,
        )

        return self.search(

            genus_taxid=genus_taxid,

            exclude=taxonomy_record.scientific_name,

        )
