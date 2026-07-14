"""
PCR product model.
"""

from dataclasses import dataclass

from envoprimer.reference.coordinates import Coordinate


@dataclass(slots=True, frozen=True)
class PCRProduct:
    """
    Represents one predicted PCR amplicon.
    """

    ####################################################################
    # Reference sequence
    ####################################################################

    chromosome: str

    ####################################################################
    # Coordinates
    ####################################################################

    forward_start: int

    reverse_end: int

    product_size: int

    forward_strand: str

    reverse_strand: str

    ####################################################################
    # Validation
    ####################################################################

    passed: bool = True

    reason: str | None = None

    ####################################################################
    # Universal coordinate
    ####################################################################

    coordinate: Coordinate | None = None

    def __post_init__(self):

        if self.coordinate is None:

            object.__setattr__(
                self,
                "coordinate",
                Coordinate(
                    start=self.forward_start,
                    end=self.reverse_end,
                    system="reference",
                ),
            )

    @property
    def length(self) -> int:

        return self.coordinate.length

    def to_dict(self) -> dict:

        return {

            "chromosome": self.chromosome,

            "forward_start": self.forward_start,

            "reverse_end": self.reverse_end,

            "product_size": self.product_size,

            "forward_strand": self.forward_strand,

            "reverse_strand": self.reverse_strand,

            "passed": self.passed,

            "reason": self.reason,

            "coordinate": self.coordinate.to_dict(),

        }

    def __repr__(self):

        return (
            f"PCRProduct("
            f"{self.chromosome}:"
            f"{self.coordinate.start}-"
            f"{self.coordinate.end}, "
            f"{self.product_size} bp)"
        )
