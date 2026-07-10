"""
Primer pair model.
"""

from dataclasses import dataclass, field

from primerforge.models.primer import Primer
from primerforge.specificity.models import SpecificityResult


@dataclass(slots=True)
class PrimerPair:
    """
    Represents a PCR primer pair.
    """

    #
    # Primer information
    #

    forward: Primer
    reverse: Primer

    product_size: int

    #
    # Overall ranking
    #

    score: float = 0.0

    breakdown: dict = field(
        default_factory=dict,
    )

    #
    # Species specificity
    #

    specificity_score: float = 100.0

    passed_specificity: bool = True

    specificity_result: SpecificityResult | None = None

    #
    # Population-aware conservation
    #

    population_conservation: float = 100.0

    population_coverage: float = 100.0

    population_result: dict = field(
        default_factory=dict,
    )
