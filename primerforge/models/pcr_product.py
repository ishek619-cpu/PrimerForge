"""
PCR product model.
"""

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class PCRProduct:
    """
    Represents one predicted PCR amplicon.
    """

    chromosome: str

    forward_start: int

    reverse_end: int

    product_size: int

    forward_strand: str

    reverse_strand: str

    passed: bool = True

    reason: str | None = None
