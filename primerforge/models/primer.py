"""
PrimerForge Primer Model
"""

from dataclasses import dataclass


@dataclass(slots=True)
class Primer:
    """
    Candidate PCR primer.
    """

    sequence: str

    start: int
    end: int

    length: int

    tm: float
    gc: float

    strand: str

    score: float = 0.0
