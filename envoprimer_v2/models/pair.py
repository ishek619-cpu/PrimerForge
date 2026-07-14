"""
Primer pair model.

Represents one candidate PCR assay.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from envoprimer_v2.models.primer import Primer


@dataclass(slots=True)
class PrimerPair:

    forward: Primer

    reverse: Primer

    product_size: int

    primer3_penalty: float = 0.0

    population_score: float = 0.0

    specificity_score: float = 0.0

    three_prime_score: float = 0.0

    thermo_score: float = 0.0

    blast_score: float = 0.0

    pcr_score: float = 0.0

    final_score: float = 0.0

    passed: bool = True

    rejection_reason: str = ""

    metadata: dict = field(default_factory=dict)

    @property
    def tm_difference(self) -> float:

        return abs(

            self.forward.tm

            - self.reverse.tm

        )

    @property
    def gc_difference(self) -> float:

        return abs(

            self.forward.gc

            - self.reverse.gc

        )

    @property
    def breakdown(self) -> dict:
        """
        Backward compatibility for legacy reports.
        """

        return {

            "primer3": self.primer3_penalty,

            "population": self.population_score,

            "specificity": self.specificity_score,

            "three_prime": self.three_prime_score,

            "thermo": self.thermo_score,

            "blast": self.blast_score,

            "pcr": self.pcr_score,

            "final": self.final_score,

        }

    def calculate_final_score(self):

        self.final_score = round(

            (

                self.population_score * 0.25 +

                self.specificity_score * 0.35 +

                self.three_prime_score * 0.20 +

                self.thermo_score * 0.10 +

                self.blast_score * 0.05 +

                self.pcr_score * 0.05

            ),

            2,

        )

        return self.final_score
