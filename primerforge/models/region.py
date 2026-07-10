"""
PrimerForge Region Model
"""

from dataclasses import dataclass


@dataclass(slots=True)
class Region:
    """
    Candidate primer design region.
    """

    snp_position: int = -1

    start: int = 0

    end: int = 0

    length: int = 0

    score: float = 0.0
