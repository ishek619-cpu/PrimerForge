"""
Closest-relative rejection engine.

Rejects primers that are predicted to amplify any
non-target species.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer_v2.specificity.scoring import (
    SpeciesSpecificityScorer,
)


@dataclass(slots=True)
class RejectionResult:

    passed: bool

    closest_species: str

    worst_identity: float

    worst_score: float

    failed_species: list[str]


class ClosestRelativeRejector:

    def __init__(self):

        self.scorer = SpeciesSpecificityScorer()

    def evaluate(

        self,

        primer: str,

        background: dict[str, str],

    ) -> RejectionResult:

        failed = []

        worst_identity = 0.0

        worst_score = 100.0

        closest = ""

        for species, sequence in background.items():

            if len(sequence) != len(primer):
                continue

            result = self.scorer.score(

                primer,

                sequence,

            )

            if result.identity > worst_identity:

                worst_identity = result.identity

                closest = species

            if result.weighted_score < worst_score:

                worst_score = result.weighted_score

            if result.passed:

                failed.append(species)

        return RejectionResult(

            passed=len(failed) == 0,

            closest_species=closest,

            worst_identity=round(
                worst_identity,
                2,
            ),

            worst_score=round(
                worst_score,
                2,
            ),

            failed_species=failed,

        )
