"""
EnvoPrimer Scoring Engine.

Combines all engine scores into one final ranking score.
"""

from __future__ import annotations


class ScoringEngine:

    def __init__(

        self,

        population_weight: float = 0.30,

        specificity_weight: float = 0.35,

        thermo_weight: float = 0.10,

        primer3_weight: float = 0.10,

        gc_weight: float = 0.05,

        pcr_weight: float = 0.10,

    ):

        self.population_weight = population_weight

        self.specificity_weight = specificity_weight

        self.thermo_weight = thermo_weight

        self.primer3_weight = primer3_weight

        self.gc_weight = gc_weight

        self.pcr_weight = pcr_weight

    def score(

        self,

        pair,

    ):

        #
        # Primer3 penalty
        #

        primer3_score = max(

            0.0,

            100.0 - pair.primer3_penalty * 10.0,

        )

        #
        # GC balance
        #

        gc_score = max(

            0.0,

            100.0 - pair.gc_difference * 5.0,

        )

        #
        # Final weighted score
        #

        final = (

            pair.population_score * self.population_weight +

            pair.specificity_score * self.specificity_weight +

            pair.thermo_score * self.thermo_weight +

            primer3_score * self.primer3_weight +

            gc_score * self.gc_weight +

            pair.pcr_score * self.pcr_weight

        )

        pair.final_score = round(

            final,

            2,

        )

        return pair.final_score
