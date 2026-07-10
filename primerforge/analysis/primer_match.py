"""
Gap-aware primer matching utilities.
"""

from __future__ import annotations

from dataclasses import dataclass

from primerforge.analysis.iupac import (
    identity,
    mismatch_count,
)


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


class PrimerMatcher:
    """
    Gap-aware primer matcher.

    Alignment coordinates are interpreted in alignment space.
    Gaps are ignored when reconstructing the biological sequence.
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
    ) -> str:
        """
        Extract primer-length sequence while ignoring gaps.
        """

        ungapped = []

        index = start

        while (
            index < len(sequence)
            and len(ungapped) < primer_length
        ):

            base = sequence[index]

            if base != "-":
                ungapped.append(base)

            index += 1

        return "".join(ungapped)

    def match(
        self,
        primer: str,
        sequence: str,
        start: int,
    ) -> PrimerMatch:

        region = self.extract_region(
            sequence,
            start,
            len(primer),
        )

        #
        # Sequence ended before primer finished.
        #
        if len(region) != len(primer):

            return PrimerMatch(
                aligned_primer=primer,
                aligned_sequence=region,
                identity=0.0,
                mismatches=len(primer),
                matched=False,
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
        )
