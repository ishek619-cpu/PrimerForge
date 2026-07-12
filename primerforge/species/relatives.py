"""
Find related species from NCBI taxonomy.
"""

from __future__ import annotations

from Bio import Entrez


class RelativeSpeciesFinder:
    """
    Find closely related species from the NCBI Taxonomy database.
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
    ) -> str:
        """
        Extract the genus from a TaxonomyRecord.
        """

        name = taxonomy_record.scientific_name

        if " " not in name:
            raise ValueError(
                f"Invalid species name: {name}"
            )

        return name.split()[0]

    def search(
        self,
        genus: str,
        exclude: str | None = None,
    ) -> list[str]:
        """
        Return all species names belonging to a genus.
        """

        handle = Entrez.esearch(
            db="taxonomy",
            term=f"{genus}[Genus] AND species[Rank]",
            retmax=1000,
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

            name = record["ScientificName"]

            if exclude is not None and name == exclude:
                continue

            species.append(name)

        species = sorted(set(species))

        return species

    def find(
        self,
        taxonomy_record,
    ) -> list[str]:
        """
        Convenience method.
        """

        genus = self.genus(
            taxonomy_record,
        )

        return self.search(
            genus=genus,
            exclude=taxonomy_record.scientific_name,
        )
