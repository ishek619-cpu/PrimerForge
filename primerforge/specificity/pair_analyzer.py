"""
Primer pair specificity analysis.
"""

from primerforge.models.pair import PrimerPair

from primerforge.specificity.validator import (
    SpecificityValidator,
)

from primerforge.specificity.pcr import (
    PCRProductFinder,
)


class PrimerPairSpecificity:
    """
    Evaluate primer pair specificity.
    """

    def __init__(self):

        self.validator = SpecificityValidator()

        self.finder = PCRProductFinder()

    def evaluate(
        self,
        pair: PrimerPair,
        database: str,
    ) -> dict:

        forward_hits = self.validator.blast_hits(
            pair.forward.sequence,
            database,
        )

        reverse_hits = self.validator.blast_hits(
            pair.reverse.sequence,
            database,
        )

        forward_score = self.validator.scorer.score(
            forward_hits,
        )

        reverse_score = self.validator.scorer.score(
            reverse_hits,
        )

        products = self.finder.find_products(
            forward_hits,
            reverse_hits,
        )

        pair_score = (
            forward_score +
            reverse_score
        ) / 2.0

        if len(products) > 1:

            pair_score -= (
                len(products) - 1
            ) * 10.0

        pair_score = round(
            max(pair_score, 0.0),
            2,
        )

        return {

            "forward_score": forward_score,

            "reverse_score": reverse_score,

            "products": products,

            "pair_score": pair_score,

            "passed": (
                pair_score >= 90.0
            ),

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

            result["pair_score"]

            / 100.0,

            2,

        )

        return pair
