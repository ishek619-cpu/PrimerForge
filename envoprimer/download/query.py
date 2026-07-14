"""
NCBI query builder for EnvoPrimer.
"""

from __future__ import annotations


MARKERS = {

    "COI": (
        "COI",
        "CO1",
        "COXI",
        "COX1",
        "\"cytochrome c oxidase subunit I\"",
        "\"cytochrome oxidase I\"",
        "\"DNA barcode\"",
        "barcode",
    ),

    "CYTB": (
        "CYTB",
        "\"cytochrome b\"",
    ),

    "12S": (
        "\"12S\"",
        "\"12S rRNA\"",
        "\"12S ribosomal RNA\"",
    ),

    "16S": (
        "\"16S\"",
        "\"16S rRNA\"",
        "\"16S ribosomal RNA\"",
    ),

    "18S": (
        "\"18S\"",
        "\"18S ribosomal RNA\"",
    ),

    "ITS": (
        "\"ITS\"",
        "\"ITS1\"",
        "\"ITS2\"",
        "\"internal transcribed spacer\"",
    ),

}


class QueryBuilder:

    def build(

        self,

        species: str,

        marker: str,

    ) -> str:

        marker = marker.upper()

        if marker not in MARKERS:

            raise ValueError(

                f"Unsupported marker: {marker}"

            )

        marker_terms = " OR ".join(

            f"{term}[All Fields]"

            for term in MARKERS[marker]

        )

        return (
            f"\"{species}\"[Organism] "
            f"AND ({marker_terms}) "
            f'NOT ("genome assembly"[Title] '
            f'OR scaffold[Title] '
            f'OR contig[Title] '
            f'OR chromosome[Title])'
        )
