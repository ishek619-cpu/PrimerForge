"""
Primer pair model.
"""

from dataclasses import dataclass, field

from primerforge.models.primer import Primer
from primerforge.specificity.models import SpecificityResult
from primerforge.thermodynamics.nearest_neighbor import (
    ThermodynamicResult,
)


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

    #
    # Nearest-neighbor thermodynamics
    #

    forward_nn: ThermodynamicResult | None = None

    reverse_nn: ThermodynamicResult | None = None

    #
    # Coverage analysis
    #

    coverage: dict = field(
        default_factory=dict,
    )

    #
    # Off-target risk
    #

    risk: dict = field(
        default_factory=dict,
    )

    #
    # Predicted PCR products
    #

    predicted_products: list = field(
        default_factory=list,
    )

    #
    # Future extensions
    #

    metadata: dict = field(
        default_factory=dict,
    )
