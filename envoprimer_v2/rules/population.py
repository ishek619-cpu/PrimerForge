"""
Population conservation rule.

Reject primer pairs that are not sufficiently conserved
across all target sequences.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult


class PopulationRule(Rule):

    name = "Population Conservation"

    mandatory = True

    def __init__(

        self,

        minimum_conservation: float = 95.0,

    ):

        self.minimum_conservation = minimum_conservation

    def evaluate(self, pair, **kwargs):

        forward = pair.forward.population_score

        reverse = pair.reverse.population_score

        conservation = min(

            forward,

            reverse,

        )

        if conservation < self.minimum_conservation:

            return RuleResult(

                name=self.name,

                passed=False,

                score=conservation,

                reason=(
                    f"Population conservation "
                    f"{conservation:.2f}% "
                    f"< {self.minimum_conservation:.2f}%"
                ),

            )

        return RuleResult(

            name=self.name,

            passed=True,

            score=round(

                conservation,

                2,

            ),

        )
