"""
EnvoPrimer Primer Model.
"""

from dataclasses import dataclass, field

from envoprimer.reference.coordinates import Coordinate


@dataclass(slots=True)
class Primer:
    """
    Represents a PCR primer.
    """

    ####################################################################
    # Sequence
    ####################################################################

    sequence: str

    ####################################################################
    # Coordinates
    ####################################################################

    start: int

    end: int

    strand: str

    ####################################################################
    # Basic properties
    ####################################################################

    length: int

    tm: float = 0.0

    gc: float = 0.0

    ####################################################################
    # Quality metrics
    ####################################################################

    gc_clamp: bool = False

    homopolymer: bool = False

    hairpin_score: float = 0.0

    self_dimer_score: float = 0.0

    score: float = 0.0

    notes: list[str] = field(
        default_factory=list,
    )

    ####################################################################
    # Universal coordinates
    ####################################################################

    coordinate: Coordinate | None = None

    ####################################################################
    # Initialization
    ####################################################################

    def __post_init__(self):

        if self.coordinate is None:

            self.coordinate = Coordinate(
                start=self.start,
                end=self.end,
                system="reference",
            )

    ####################################################################
    # Convenience properties
    ####################################################################

    @property
    def coordinate_length(self) -> int:

        return self.coordinate.length

    def to_dict(self) -> dict:

        return {

            "sequence": self.sequence,

            "start": self.start,

            "end": self.end,

            "strand": self.strand,

            "length": self.length,

            "tm": self.tm,

            "gc": self.gc,

            "coordinate": self.coordinate.to_dict(),

        }

    def __repr__(self):

        return (
            f"Primer("
            f"{self.sequence}, "
            f"{self.coordinate})"
        )
