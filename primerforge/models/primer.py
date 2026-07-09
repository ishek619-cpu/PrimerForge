"""
PrimerForge Primer Model
"""

from dataclasses import dataclass, field


@dataclass(slots=True)
class Primer:
    """
    Represents a PCR primer.
    """

    sequence: str

    start: int
    end: int

    strand: str

    length: int

    tm: float = 0.0
    gc: float = 0.0

    gc_clamp: bool = False

    homopolymer: bool = False

    hairpin_score: float = 0.0

    self_dimer_score: float = 0.0

    score: float = 0.0

    notes: list[str] = field(default_factory=list)
