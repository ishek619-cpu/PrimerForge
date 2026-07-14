"""
Diagnostic SNP rule.

Requires at least one primer to terminate on, or very near,
a diagnostic SNP.
"""

from __future__ import annotations

from envoprimer_v2.rules.base import Rule, RuleResult
from envoprimer_v2.diagnostics.three_prime import ThreePrimeAnalyzer


class DiagnosticSNPRule(Rule):

    name = "Diagnostic SNP"

    mandatory = True

    def __init__(

        self,

        maximum_distance: int = 2,

    ):

        self.analyzer = ThreePrimeAnalyzer(

            maximum_distance=maximum_distance,

        )

    def evaluate(

        self,

        pair,

        diagnostic_snps,

    ):

        forward = self.analyzer.evaluate(

            primer_start=pair.forward.start,

            primer_end=pair.forward.end,

            strand="+",

            diagnostic_snps=diagnostic_snps,

        )

        reverse = self.analyzer.evaluate(

            primer_start=pair.reverse.start,

            primer_end=pair.reverse.end,

            strand="-",

            diagnostic_snps=diagnostic_snps,

        )

        if not (forward.passed or reverse.passed):

            return RuleResult(

                name=self.name,

                passed=False,

                score=0.0,

                reason="No diagnostic SNP within 2 bp of either primer 3' end.",

            )

        score = max(

            forward.score,

            reverse.score,

        )

        return RuleResult(

            name=self.name,

            passed=True,

            score=round(score, 2),

        )
