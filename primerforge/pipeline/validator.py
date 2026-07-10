"""
Primer validation pipeline.

Coordinates all validation modules.
"""

from pathlib import Path

from primerforge.coverage.analyzer import CoverageAnalyzer
from primerforge.pcr.product_predictor import PCRProductPredictor
from primerforge.primer.pairscore import PrimerPairScorer
from primerforge.risk.assessor import RiskAssessor
from primerforge.specificity.engine import SpecificityEngine
from primerforge.specificity.three_prime import ThreePrimeAnalyzer
from primerforge.thermodynamics.nearest_neighbor import (
    NearestNeighborCalculator,
)


class PrimerValidationPipeline:
    """
    Complete validation workflow.
    """

    def __init__(self):

        self.coverage = CoverageAnalyzer()

        self.predictor = PCRProductPredictor()

        self.risk = RiskAssessor()

        self.three_prime = ThreePrimeAnalyzer()

        self.scorer = PrimerPairScorer()

        self.nn = NearestNeighborCalculator()

    def validate(
        self,
        pairs,
        reference_fasta: Path,
        alignment_fasta: Path | None = None,
        specificity_engine: SpecificityEngine | None = None,
    ):

        validated = []

        for pair in pairs:

            #
            # Thermodynamics
            #
            pair.forward_nn = self.nn.calculate(
                pair.forward.sequence,
            )

            pair.reverse_nn = self.nn.calculate(
                pair.reverse.sequence,
            )

            #
            # Coverage
            #
            if alignment_fasta is not None:

                pair.coverage = self.coverage.analyse(
                    pair,
                    alignment_fasta,
                )

            #
            # Species specificity
            #
            if specificity_engine is not None:

                result = specificity_engine.evaluate_pair(
                    pair,
                    Path("results/blast"),
                )

                for hit in result.off_target_hits:

                    hit.three_prime_mismatches = (
                        self.three_prime.count(hit)
                    )

                    hit.three_prime_score = (
                        self.three_prime.score(hit)
                    )

                    risk = self.risk.assess(hit)

                    if (
                        not hasattr(pair, "risk")
                        or risk["score"] > pair.risk["score"]
                    ):
                        pair.risk = risk

            #
            # Initial scoring
            #
            validation = {
                "thermo": 100.0,
                "conservation": 100.0,
                "snp": 100.0,
            }

            self.scorer.score(
                pair,
                validation,
            )

            validated.append(pair)

        return self.scorer.rank(validated)
