"""
Data models for species-specific primer analysis.
"""

from dataclasses import dataclass, field
from typing import List


@dataclass(frozen=True)
class BlastHit:
    """
    Single BLAST alignment hit.
    """

    accession: str

    species: str

    identity: float

    coverage: float

    alignment_length: int

    mismatches: int

    gap_opens: int

    qstart: int

    qend: int

    sstart: int

    send: int

    #
    # NEW
    #
    query_sequence: str = ""

    subject_sequence: str = ""

    strand: str = "plus"

    bitscore: float = 0.0

    evalue: float = 1.0

    #
    # Filled later by ThreePrimeAnalyzer
    #
    three_prime_mismatches: int = 0

    three_prime_score: float = 100.0


@dataclass(frozen=True)
class OffTargetHit(BlastHit):
    """
    Off-target BLAST hit.
    """

    penalty: float = 0.0


@dataclass(frozen=True)
class SpecificityResult:
    """
    Final specificity evaluation.
    """

    forward_hits: List[BlastHit] = field(
        default_factory=list,
    )

    reverse_hits: List[BlastHit] = field(
        default_factory=list,
    )

    target_hits: List[BlastHit] = field(
        default_factory=list,
    )

    off_target_hits: List[OffTargetHit] = field(
        default_factory=list,
    )

    specificity_score: float = 100.0

    passed: bool = True

    rejection_reason: str | None = None
