"""
Sequence metadata used throughout EnvoPrimer.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class SequenceMetadata:
    """
    Metadata describing one downloaded sequence.
    """

    accession: str

    species: str

    marker: str

    definition: str

    length: int

    source: str = "NCBI"

    role: str = ""

    fasta: str = ""
