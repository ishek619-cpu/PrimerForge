"""
Homopolymer rule.

Reject primers containing long homopolymer runs.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult


class HomopolymerRule(Rule):

    name = "Homopolymer"

    mandatory = True

    def __init__(

        self,

        maximum_run: int = 4,

    ):

        self.maximum_run = maximum_run

    def longest_run(self, sequence: str) -> int:

        if not sequence:
            return 0

        longest = 1
        current = 1

        for previous, current_base in zip(sequence, sequence[1:]):

            if previous == current_base:

                current += 1

                if current > longest:

                    longest = current

            else:

                current = 1

        return longest

    def evaluate(self, pair, **kwargs):

        forward_run = self.longest_run(

            pair.forward.sequence.upper(),

        )

        reverse_run = self.longest_run(

            pair.reverse.sequence.upper(),

        )

        longest = max(

            forward_run,

            reverse_run,

        )

        if longest > self.maximum_run:

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason=f"Homopolymer length {longest} exceeds {self.maximum_run}.",

            )

        score = max(

            0.0,

            100.0 - (longest - 1) * 10.0,

        )

        return RuleResult(

            name=self.name,

            passed=True,

            score=round(score, 2),

        )
