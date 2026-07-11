"""
Gap-aware primer matching utilities.
"""

from __future__ import annotations

from dataclasses import dataclass

from primerforge.analysis.iupac import (
    identity,
    mismatch_count,
)
from primerforge.reference.coordinates import Coordinate


@dataclass(frozen=True, slots=True)
class PrimerMatch:
    """
    Result of matching one primer against one sequence.
    """

    aligned_primer: str

    aligned_sequence: str

    identity: float

    mismatches: int

    matched: bool

    coordinate: Coordinate | None = None


class PrimerMatcher:
    """
    Gap-aware primer matcher.

    Coordinates are interpreted in alignment space.
    """

    def __init__(
        self,
        max_mismatches: int = 1,
    ):

        self.max_mismatches = max_mismatches

    def extract_region(
        self,
        sequence: str,
        start: int,
        primer_length: int,
    ) -> tuple[str, Coordinate]:

        ungapped = []

        index = start

        alignment_start = index + 1

        alignment_end = alignment_start

        while (
            index < len(sequence)
            and len(ungapped) < primer_length
        ):

            base = sequence[index]

            if base != "-":

                ungapped.append(
                    base,
                )

                alignment_end = index + 1

            index += 1

        coordinate = Coordinate(

            start=alignment_start,

            end=alignment_end,

            system="alignment",

        )

        return "".join(
            ungapped,
        ), coordinate

    def match(
        self,
        primer: str,
        sequence: str,
        start: int,
    ) -> PrimerMatch:

        region, coordinate = self.extract_region(

            sequence,

            start,

            len(primer),

        )

        if len(region) != len(primer):

            return PrimerMatch(

                aligned_primer=primer,

                aligned_sequence=region,

                identity=0.0,

                mismatches=len(primer),

                matched=False,

                coordinate=coordinate,

            )

        mm = mismatch_count(

            primer,

            region,

        )

        ident = identity(

            primer,

            region,

        )

        return PrimerMatch(

            aligned_primer=primer,

            aligned_sequence=region,

            identity=ident,

            mismatches=mm,

            matched=mm <= self.max_mismatches,

            coordinate=coordinate,

        )
