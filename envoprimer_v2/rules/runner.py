"""
EnvoPrimer Rule Engine.

Evaluates every primer pair using the complete
species-specific rule set.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from envoprimer_v2.rules.base import RuleResult


@dataclass(slots=True)
class RuleEngineResult:

    passed: bool

    total_score: float

    failed_rule: str = ""

    failed_reason: str = ""

    results: list[RuleResult] = field(default_factory=list)


class RuleEngine:

    def __init__(

        self,

        debug: bool = False,

    ):

        self.rules = []

        self.debug = debug

    def add(self, rule):

        self.rules.append(rule)

    def evaluate(

        self,

        pair,

        **context,

    ) -> RuleEngineResult:

        results = []

        scores = []

        failed_rule = ""

        failed_reason = ""

        passed = True

        if self.debug:

            print()
            print("=" * 80)
            print("Evaluating Primer Pair")
            print("=" * 80)
            print("Forward :", pair.forward.sequence)
            print("Reverse :", pair.reverse.sequence)
            print("Product :", pair.product_size)
            print()

        for rule in self.rules:

            result = rule.evaluate(

                pair,

                **context,

            )

            results.append(result)

            scores.append(result.score)

            if self.debug:

                state = "PASS" if result.passed else "FAIL"

                print(

                    f"{rule.name:<30}"

                    f"{state:<6}"

                    f"{result.score:>7.2f}"

                )

                if result.reason:

                    print("   ", result.reason)

            if rule.mandatory and not result.passed:

                passed = False

                failed_rule = rule.name

                failed_reason = result.reason

                break

        total = 0.0

        if passed and scores:

            total = round(

                sum(scores)

                / len(scores),

                2,

            )

        if self.debug:

            print()

            if passed:

                print(

                    f"ACCEPTED (Score={total})"

                )

            else:

                print(

                    f"REJECTED ({failed_rule})"

                )

            print("=" * 80)

        return RuleEngineResult(

            passed=passed,

            total_score=total,

            failed_rule=failed_rule,

            failed_reason=failed_reason,

            results=results,

        )
