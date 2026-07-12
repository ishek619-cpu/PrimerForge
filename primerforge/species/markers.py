"""
Supported genetic markers and their GenBank search queries.
"""

from __future__ import annotations

MARKERS = {

    "COI": [
        "COI",
        "CO1",
        "COX1",
        "COXI",
        "cytochrome oxidase I",
        "cytochrome c oxidase subunit I",
        "DNA barcode",
        "barcode",
    ],

    "CYTB": [
        "CYTB",
        "cytb",
        "cytochrome b",
        "CYTOCHROME B",
    ],

    "12S": [
        "12S",
        "12S rRNA",
        "12S ribosomal RNA",
    ],

    "16S": [
        "16S",
        "16S rRNA",
        "16S ribosomal RNA",
    ],

    "18S": [
        "18S",
        "18S rRNA",
        "18S ribosomal RNA",
    ],

    "ITS": [
        "ITS",
        "ITS1",
        "ITS2",
        "internal transcribed spacer",
    ],
}


def marker_synonyms(marker: str) -> list[str]:

    marker = marker.upper()

    if marker not in MARKERS:

        raise KeyError(
            f"Unsupported marker: {marker}"
        )

    return MARKERS[marker]


def build_query(
    species: str,
    marker: str,
) -> str:

    synonyms = marker_synonyms(
        marker,
    )

    marker_query = " OR ".join(
        f'"{name}"[Title]'
        for name in synonyms
    )

    return (
        f'"{species}"[Organism] '
        f'AND ({marker_query}) '
        f'AND mitochondrion[Filter] '
        f'NOT ("whole genome"[Title] '
        f'OR chromosome[Title] '
        f'OR scaffold[Title] '
        f'OR contig[Title] '
        f'OR WGS[Title])'
    )
