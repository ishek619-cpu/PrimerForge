"""
3' end quality rule.

Rejects primers with weak or unstable 3' ends.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult


class ThreePrimeRule(Rule):

    name = "3' End"

    mandatory = True

    def __init__(

        self,

        tail_length: int = 5,

        minimum_gc: int = 1,

        maximum_at: int = 4,

    ):

        self.tail_length = tail_length
        self.minimum_gc = minimum_gc
        self.maximum_at = maximum_at

    def evaluate_primer(self, sequence: str):

        tail = sequence.upper()[-self.tail_length:]

        gc = sum(

            base in {"G", "C"}

            for base in tail

        )

        at = self.tail_length - gc

        if gc < self.minimum_gc:

            return False, 0.0, "Insufficient GC at 3' end."

        if at > self.maximum_at:

            return False, 0.0, "AT-rich 3' end."

        score = 100.0

        score -= abs(gc - 2) * 10.0

        return True, max(score, 0.0), ""

    def evaluate(self, pair, **kwargs):

        ok_f, score_f, reason_f = self.evaluate_primer(

            pair.forward.sequence,

        )

        if not ok_f:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason="Forward: " + reason_f,

            )

        ok_r, score_r, reason_r = self.evaluate_primer(

            pair.reverse.sequence,

        )

        if not ok_r:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason="Reverse: " + reason_r,

            )

        return RuleResult(

            name=self.name,

            passed=True,

            score=round(

                (score_f + score_r) / 2,

                2,

            ),

        )
