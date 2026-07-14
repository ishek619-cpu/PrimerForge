"""
Primer validation pipeline.

Coordinates all validation modules.
"""

from pathlib import Path

from envoprimer.coverage.analyzer import CoverageAnalyzer
from envoprimer.pcr.product_predictor import PCRProductPredictor
from envoprimer.primer.pairscore import PrimerPairScorer
from envoprimer.risk.assessor import RiskAssessor
from envoprimer.specificity.engine import SpecificityEngine
from envoprimer.specificity.three_prime import ThreePrimeAnalyzer
from envoprimer.thermodynamics.nearest_neighbor import (
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

                #
                # Assess risk without modifying frozen OffTargetHit objects
                #
                highest_risk = None

                for hit in result.off_target_hits:

                    risk = self.risk.assess(hit)

                    if (
                        highest_risk is None
                        or risk["score"] > highest_risk["score"]
                    ):
                        highest_risk = risk

                if highest_risk is not None:
                    pair.risk = highest_risk

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
