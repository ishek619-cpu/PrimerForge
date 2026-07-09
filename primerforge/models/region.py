"""
PrimerForge Region Model
"""

from dataclasses import dataclass


@dataclass(slots=True)
class Region:
    """
    Candidate primer design region.
    """

    snp_position: int
    start: int
    end: int
    length: int
