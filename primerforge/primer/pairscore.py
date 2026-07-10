"""
Primer pair scoring engine.
"""

from primerforge.models.pair import PrimerPair
from primerforge.primer3.thermo import ThermoAnalyzer
from primerforge.specificity.validator import SpecificityValidator
from primerforge.specificity.pair_analyzer import (
    PrimerPairSpecificity,
)


class PrimerPairScorer:
    """
    Score and rank primer pairs.
    """

    def __init__(self):

        self.thermo = ThermoAnalyzer()

        self.validator = SpecificityValidator()

        self.pair_specificity = PrimerPairSpecificity()

    def score(
        self,
        pair: PrimerPair,
        database: str | None = None,
    ) -> PrimerPair:

        pair.forward = self.thermo.evaluate(
            pair.forward,
        )

        pair.reverse = self.thermo.evaluate(
            pair.reverse,
        )

        heterodimer = self.thermo.heterodimer(
            pair.forward,
            pair.reverse,
        )

        score = 100.0

        dtm = abs(
            pair.forward.tm -
            pair.reverse.tm
        )

        score -= dtm * 5.0

        dgc = abs(
            pair.forward.gc -
            pair.reverse.gc
        )

        score -= dgc * 0.5

        if not (
            80 <= pair.product_size <= 250
        ):
            score -= 20

        if pair.forward.gc_clamp:
            score += 2

        if pair.reverse.gc_clamp:
            score += 2

        score -= pair.forward.hairpin_score
        score -= pair.reverse.hairpin_score

        score -= pair.forward.self_dimer_score
        score -= pair.reverse.self_dimer_score

        score -= heterodimer

        if database is not None:

            pair.score = round(
                max(score, 0.0),
                2,
            )

            pair = self.pair_specificity.annotate(
                pair,
                database,
            )

        else:

            pair.score = round(
                max(score, 0.0),
                2,
            )

        return pair

    def rank(
        self,
        pairs,
        database: str | None = None,
    ):

        scored = [

            self.score(
                pair,
                database,
            )

            for pair in pairs

        ]

        return sorted(
            scored,
            key=lambda pair: pair.score,
            reverse=True,
        )
