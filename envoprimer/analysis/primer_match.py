"""
Gap-aware primer matching utilities.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer.analysis.iupac import compare
from envoprimer.reference.coordinates import Coordinate


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

    Coordinates supplied to this class are 1-based alignment
    coordinates. Python strings are 0-based, so we convert them.
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

        # Convert 1-based alignment coordinate to Python index
        index = start - 1

        alignment_start = start
        alignment_end = start

        ungapped = []

        while (
            index < len(sequence)
            and len(ungapped) < primer_length
        ):

            base = sequence[index]

            if base != "-":
                ungapped.append(base)
                alignment_end = index + 1

            index += 1

        coordinate = Coordinate(
            start=alignment_start,
            end=alignment_end,
            system="alignment",
        )

        return "".join(ungapped), coordinate

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

        mismatches, identity = compare(
            primer,
            region,
        )

        return PrimerMatch(
            aligned_primer=primer,
            aligned_sequence=region,
            identity=identity,
            mismatches=mismatches,
            matched=mismatches <= self.max_mismatches,
            coordinate=coordinate,
        )
