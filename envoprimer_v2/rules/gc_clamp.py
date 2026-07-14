"""
GC clamp rule.

Ensures a stable 3' GC clamp without making it excessively GC-rich.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult


class GCClampRule(Rule):

    name = "GC Clamp"

    mandatory = True

    def __init__(

        self,

        minimum_gc_bases: int = 1,

        maximum_gc_bases: int = 3,

    ):

        self.minimum_gc_bases = minimum_gc_bases
        self.maximum_gc_bases = maximum_gc_bases

    def gc_clamp(self, sequence: str) -> int:

        tail = sequence.upper()[-5:]

        return sum(

            base in {"G", "C"}

            for base in tail

        )

    def evaluate(self, pair, **kwargs):

        forward_gc = self.gc_clamp(

            pair.forward.sequence,

        )

        reverse_gc = self.gc_clamp(

            pair.reverse.sequence,

        )

        if (

            forward_gc < self.minimum_gc_bases

            or

            forward_gc > self.maximum_gc_bases

        ):

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason=f"Forward primer GC clamp ({forward_gc}/5) outside allowed range.",

            )

        if (

            reverse_gc < self.minimum_gc_bases

            or

            reverse_gc > self.maximum_gc_bases

        ):

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason=f"Reverse primer GC clamp ({reverse_gc}/5) outside allowed range.",

            )

        score = 100.0

        if forward_gc == 2:

            score += 2

        if reverse_gc == 2:

            score += 2

        return RuleResult(

            name=self.name,

            passed=True,

            score=min(100.0, score),

        )
