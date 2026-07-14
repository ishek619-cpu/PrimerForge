"""
Primer pair ranking.
"""

from envoprimer.models.pair import PrimerPair
from envoprimer.thermodynamics.nearest_neighbor import (
    NearestNeighborCalculator,
)


class PrimerPairScorer:
    """
    Multi-criteria primer pair scoring.
    """

    def __init__(self):

        self.nn = NearestNeighborCalculator()

    def score(
        self,
        pair: PrimerPair,
        validation: dict,
        specificity: dict | None = None,
        multiplex: dict | None = None,
    ) -> float:

        thermo = validation.get(
            "thermo",
            100.0,
        )

        conservation = validation.get(
            "conservation",
            100.0,
        )

        snp = validation.get(
            "snp",
            100.0,
        )

        #
        # Product size
        #
        if 80 <= pair.product_size <= 150:
            product = 100.0
        elif 151 <= pair.product_size <= 200:
            product = 90.0
        elif 60 <= pair.product_size <= 250:
            product = 80.0
        else:
            product = 50.0

        #
        # Tm balance
        #
        tm_difference = abs(
            pair.forward.tm
            - pair.reverse.tm
        )

        tm_balance = max(
            0.0,
            100.0 - tm_difference * 20.0,
        )

        #
        # GC balance
        #
        gc_difference = abs(
            pair.forward.gc
            - pair.reverse.gc
        )

        gc_balance = max(
            0.0,
            100.0 - gc_difference * 5.0,
        )

        #
        # Species specificity
        #
        if pair.passed_specificity:
            specificity_score = pair.specificity_score
        else:
            specificity_score = 0.0

        if specificity is not None:
            specificity_score = specificity.get(
                "specificity",
                specificity_score,
            )

        #
        # Multiplex
        #
        multiplex_score = (
            multiplex.get(
                "compatibility",
                100.0,
            )
            if multiplex
            else 100.0
        )

        #
        # Nearest-neighbor thermodynamics
        #
        forward_nn = self.nn.calculate(
            pair.forward.sequence,
        )

        reverse_nn = self.nn.calculate(
            pair.reverse.sequence,
        )

        tm_nn_difference = abs(
            forward_nn.tm
            - reverse_nn.tm
        )

        nn_score = max(
            0.0,
            100.0 - tm_nn_difference,
        )

        #
        # Final weighted score
        #
        final_score = (

            thermo * 0.20 +

            conservation * 0.20 +

            snp * 0.15 +

            specificity_score * 0.15 +

            nn_score * 0.10 +

            product * 0.10 +

            tm_balance * 0.05 +

            gc_balance * 0.03 +

            multiplex_score * 0.02

        )

        pair.score = round(
            final_score,
            2,
        )

        pair.breakdown = {

            "thermo": round(
                thermo,
                2,
            ),

            "conservation": round(
                conservation,
                2,
            ),

            "snp": round(
                snp,
                2,
            ),

            "specificity": round(
                specificity_score,
                2,
            ),

            "nearest_neighbor": round(
                nn_score,
                2,
            ),

            "product": round(
                product,
                2,
            ),

            "tm_balance": round(
                tm_balance,
                2,
            ),

            "gc_balance": round(
                gc_balance,
                2,
            ),

            "multiplex": round(
                multiplex_score,
                2,
            ),

            "final_score": pair.score,

        }

        return pair.score

    def rank(
        self,
        pairs: list[PrimerPair],
    ) -> list[PrimerPair]:

        return sorted(
            pairs,
            key=lambda pair: pair.score,
            reverse=True,
        )
