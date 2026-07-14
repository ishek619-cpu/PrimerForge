"""
Melting temperature rule.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult


class TmRule(Rule):

    name = "Melting Temperature"

    mandatory = True

    def __init__(

        self,

        minimum: float = 58.0,

        optimum: float = 60.0,

        maximum: float = 62.0,

        maximum_difference: float = 2.0,

    ):

        self.minimum = minimum
        self.optimum = optimum
        self.maximum = maximum
        self.maximum_difference = maximum_difference

    def evaluate(self, pair, **kwargs):

        forward_tm = pair.forward.tm
        reverse_tm = pair.reverse.tm

        if forward_tm < self.minimum or forward_tm > self.maximum:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason="Forward primer Tm outside acceptable range.",

            )

        if reverse_tm < self.minimum or reverse_tm > self.maximum:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason="Reverse primer Tm outside acceptable range.",

            )

        difference = abs(

            forward_tm

            - reverse_tm

        )

        if difference > self.maximum_difference:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason=f"Tm difference {difference:.2f}°C exceeds {self.maximum_difference:.2f}°C.",

            )

        mean_tm = (

            forward_tm

            + reverse_tm

        ) / 2.0

        score = max(

            0.0,

            100.0

            - abs(mean_tm - self.optimum) * 15.0

            - difference * 10.0,

        )

        return RuleResult(

            name=self.name,

            passed=True,

            score=round(score, 2),

        )
