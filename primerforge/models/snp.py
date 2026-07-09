"""
PrimerForge SNP Model
"""

from dataclasses import dataclass


@dataclass(slots=True)
class SNP:
    """
    Represents a single SNP in an alignment.
    """

    position: int
    reference: str
    alternatives: list[str]
    counts: dict[str, int]
