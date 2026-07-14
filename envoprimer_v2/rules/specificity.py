"""
Species specificity rule.

Rejects primer pairs predicted to amplify the closest
non-target species.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult
from envoprimer_v2.specificity.rejection import (
    ClosestRelativeRejector,
)


class SpeciesSpecificityRule(Rule):

    name = "Species Specificity"

    mandatory = True

    def __init__(self):

        self.rejector = ClosestRelativeRejector()

    def evaluate(

        self,

        pair,

        background_sequences: dict[str, tuple[str, str]],

    ):

        forward_background = {

            species: sequences[0]

            for species, sequences in background_sequences.items()

        }

        reverse_background = {

            species: sequences[1]

            for species, sequences in background_sequences.items()

        }

        forward = self.rejector.evaluate(

            pair.forward.sequence,

            forward_background,

        )

        reverse = self.rejector.evaluate(

            pair.reverse.sequence,

            reverse_background,

        )

        if not forward.passed:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason=(
                    "Forward primer amplifies "
                    f"{forward.closest_species}"
                ),

            )

        if not reverse.passed:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason=(
                    "Reverse primer amplifies "
                    f"{reverse.closest_species}"
                ),

            )

        score = min(

            forward.worst_score,

            reverse.worst_score,

        )

        return RuleResult(

            name=self.name,

            passed=True,

            score=round(score, 2),

        )
