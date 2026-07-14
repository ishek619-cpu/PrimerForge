"""
Species-specific primer scoring.

Weighted mismatch scoring with strong emphasis on the 3' end.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class SpecificityScore:

    identity: float

    mismatches: int

    weighted_score: float

    three_prime_mismatches: int

    passed: bool


class SpeciesSpecificityScorer:

    def __init__(

        self,

        maximum_identity: float = 85.0,

    ):

        self.maximum_identity = maximum_identity

    def score(

        self,

        primer: str,

        target: str,

    ) -> SpecificityScore:

        primer = primer.upper()

        target = target.upper()

        if len(primer) != len(target):

            raise ValueError(
                "Primer and target must have equal length."
            )

        matches = 0

        mismatches = 0

        weighted = 0.0

        three_prime = 0

        length = len(primer)

        for i, (a, b) in enumerate(zip(primer, target)):

            if a == b:

                matches += 1
                continue

            mismatches += 1

            distance = length - i

            #
            # Strong weighting near the 3' end
            #
            if distance <= 5:

                weighted += 3.0
                three_prime += 1

            elif distance <= 10:

                weighted += 2.0

            else:

                weighted += 1.0

        identity = matches / length * 100.0

        specificity = 100.0 - (

            weighted

            / (length * 3.0)

            * 100.0

        )

        passed = (

            identity <= self.maximum_identity

            or

            three_prime > 0

        )

        return SpecificityScore(

            identity=round(identity, 2),

            mismatches=mismatches,

            weighted_score=round(

                specificity,

                2,

            ),

            three_prime_mismatches=three_prime,

            passed=passed,

        )
