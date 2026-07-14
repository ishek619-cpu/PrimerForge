"""
Primer pair model.
"""

from dataclasses import dataclass, field

from envoprimer.models.primer import Primer
from envoprimer.reference.coordinates import Coordinate
from envoprimer.specificity.models import SpecificityResult
from envoprimer.thermodynamics.nearest_neighbor import (
    ThermodynamicResult,
)


@dataclass(slots=True)
class PrimerPair:
    """
    Represents a PCR primer pair.
    """

    ####################################################################
    # Primers
    ####################################################################

    forward: Primer

    reverse: Primer

    product_size: int

    ####################################################################
    # Universal coordinate
    ####################################################################

    coordinate: Coordinate | None = None

    ####################################################################
    # Overall ranking
    ####################################################################

    score: float = 0.0

    breakdown: dict = field(
        default_factory=dict,
    )

    ####################################################################
    # Species specificity
    ####################################################################

    specificity_score: float = 100.0

    passed_specificity: bool = True

    specificity_result: SpecificityResult | None = None

    ####################################################################
    # Population conservation
    ####################################################################

    population_conservation: float = 100.0

    population_coverage: float = 100.0

    population_result: dict = field(
        default_factory=dict,
    )

    ####################################################################
    # Thermodynamics
    ####################################################################

    forward_nn: ThermodynamicResult | None = None

    reverse_nn: ThermodynamicResult | None = None

    ####################################################################
    # Coverage
    ####################################################################

    coverage: dict = field(
        default_factory=dict,
    )

    ####################################################################
    # Risk assessment
    ####################################################################

    risk: dict = field(
        default_factory=dict,
    )

    ####################################################################
    # PCR products
    ####################################################################

    predicted_products: list = field(
        default_factory=list,
    )

    ####################################################################
    # Metadata
    ####################################################################

    metadata: dict = field(
        default_factory=dict,
    )

    ####################################################################
    # Initialization
    ####################################################################

    def __post_init__(self):

        if self.coordinate is None:

            self.coordinate = Coordinate(

                start=self.forward.coordinate.start,

                end=self.reverse.coordinate.end,

                system=self.forward.coordinate.system,

            )

    ####################################################################
    # Convenience
    ####################################################################

    @property
    def length(self) -> int:

        return self.coordinate.length

    def to_dict(self) -> dict:

        return {

            "forward": self.forward.to_dict(),

            "reverse": self.reverse.to_dict(),

            "product_size": self.product_size,

            "score": self.score,

            "population_conservation": self.population_conservation,

            "population_coverage": self.population_coverage,

            "specificity_score": self.specificity_score,

            "coordinate": self.coordinate.to_dict(),

        }

    def __repr__(self):

        return (
            f"PrimerPair("
            f"{self.forward.sequence} / "
            f"{self.reverse.sequence}, "
            f"{self.coordinate.start}-"
            f"{self.coordinate.end})"
        )
