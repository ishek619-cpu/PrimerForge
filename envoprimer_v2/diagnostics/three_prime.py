"""
3' diagnostic SNP enforcement.

A primer should terminate on, or immediately adjacent to,
a diagnostic SNP whenever possible. This greatly improves
species specificity by destabilizing extension in non-target
species.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer_v2.diagnostics.snp import DiagnosticSNP


@dataclass(slots=True)
class ThreePrimeResult:

    passed: bool

    distance: int

    score: float

    snp: DiagnosticSNP | None


class ThreePrimeAnalyzer:

    def __init__(

        self,

        maximum_distance: int = 2,

    ):

        self.maximum_distance = maximum_distance

    def evaluate(

        self,

        primer_start: int,

        primer_end: int,

        strand: str,

        diagnostic_snps: list[DiagnosticSNP],

    ) -> ThreePrimeResult:

        if not diagnostic_snps:

            return ThreePrimeResult(

                passed=False,

                distance=-1,

                score=0.0,

                snp=None,

            )

        if strand == "+":

            three_prime = primer_end

        else:

            three_prime = primer_start

        nearest = None

        best_distance = 10**9

        for snp in diagnostic_snps:

            distance = abs(

                snp.alignment_position

                - three_prime

            )

            if distance < best_distance:

                best_distance = distance

                nearest = snp

        passed = (

            best_distance

            <= self.maximum_distance

        )

        if passed:

            score = max(

                0.0,

                100.0 - best_distance * 20.0,

            )

        else:

            score = 0.0

        return ThreePrimeResult(

            passed=passed,

            distance=best_distance,

            score=round(score, 2),

            snp=nearest,

        )
