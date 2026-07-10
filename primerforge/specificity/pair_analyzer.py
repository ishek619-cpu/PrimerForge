"""
Primer pair specificity analysis.
"""

from primerforge.models.pair import PrimerPair
from primerforge.specificity.validator import (
    SpecificityValidator,
)


class PrimerPairSpecificity:
    """
    Evaluate primer pair specificity.
    """

    def __init__(self):

        self.validator = SpecificityValidator()

    def evaluate(
        self,
        pair: PrimerPair,
        database: str,
    ) -> dict:

        forward_score = self.validator.validate(
            pair.forward.sequence,
            database,
        )

        reverse_score = self.validator.validate(
            pair.reverse.sequence,
            database,
        )

        pair_score = round(
            (
                forward_score +
                reverse_score
            ) / 2.0,
            2,
        )

        return {

            "forward_score": forward_score,

            "reverse_score": reverse_score,

            "pair_score": pair_score,

            "passed": pair_score >= 90.0,

        }

    def annotate(
        self,
        pair: PrimerPair,
        database: str,
    ) -> PrimerPair:

        result = self.evaluate(
            pair,
            database,
        )

        pair.score = round(
            pair.score *
            result["pair_score"] / 100.0,
            2,
        )

        return pair
