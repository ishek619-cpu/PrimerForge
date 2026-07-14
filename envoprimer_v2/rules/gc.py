"""
GC content rule.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult


class GCRule(Rule):

    name = "GC Content"

    mandatory = True

    def __init__(

        self,

        minimum: float = 40.0,

        optimum: float = 50.0,

        maximum: float = 60.0,

        maximum_difference: float = 10.0,

    ):

        self.minimum = minimum
        self.optimum = optimum
        self.maximum = maximum
        self.maximum_difference = maximum_difference

    def evaluate(self, pair, **kwargs):

        forward_gc = pair.forward.gc
        reverse_gc = pair.reverse.gc

        if forward_gc < self.minimum or forward_gc > self.maximum:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason="Forward primer GC outside acceptable range.",

            )

        if reverse_gc < self.minimum or reverse_gc > self.maximum:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason="Reverse primer GC outside acceptable range.",

            )

        difference = abs(

            forward_gc

            - reverse_gc

        )

        if difference > self.maximum_difference:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason=f"GC difference {difference:.2f}% exceeds {self.maximum_difference:.2f}%",

            )

        mean_gc = (

            forward_gc

            + reverse_gc

        ) / 2.0

        score = max(

            0.0,

            100.0

            - abs(mean_gc - self.optimum) * 3.0

            - difference * 2.0,

        )

        return RuleResult(

            name=self.name,

            passed=True,

            score=round(score, 2),

        )
