"""
Primer pair model.
"""

from dataclasses import dataclass

from primerforge.models.primer import Primer


@dataclass(slots=True)
class PrimerPair:
    """
    Represents a PCR primer pair.
    """

    forward: Primer
    reverse: Primer

    product_size: int

    score: float = 0.0
