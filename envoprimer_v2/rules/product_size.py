"""
Product size rule.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult


class ProductSizeRule(Rule):

    name = "Product Size"

    mandatory = True

    def __init__(
        self,
        minimum: int = 70,
        optimum_low: int = 120,
        optimum_high: int = 150,
        maximum: int = 200,
    ):

        self.minimum = minimum
        self.optimum_low = optimum_low
        self.optimum_high = optimum_high
        self.maximum = maximum

    def evaluate(
        self,
        pair,
        **kwargs,
    ):

        size = pair.product_size

        if size < self.minimum:

            return RuleResult(
                name=self.name,
                passed=False,
                score=0.0,
                reason=f"Amplicon {size} bp < {self.minimum} bp",
            )

        if size > self.maximum:

            return RuleResult(
                name=self.name,
                passed=False,
                score=0.0,
                reason=f"Amplicon {size} bp > {self.maximum} bp",
            )

        if self.optimum_low <= size <= self.optimum_high:

            score = 100.0

        elif size < self.optimum_low:

            score = max(
                60.0,
                100.0 - (self.optimum_low - size) * 0.75,
            )

        else:

            score = max(
                60.0,
                100.0 - (size - self.optimum_high) * 0.75,
            )

        return RuleResult(
            name=self.name,
            passed=True,
            score=round(score, 2),
        )
