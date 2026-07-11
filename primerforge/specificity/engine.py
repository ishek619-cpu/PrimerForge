"""
Species-specific primer evaluation engine.
"""

from pathlib import Path

from primerforge.models.pair import PrimerPair
from primerforge.reference.coordinates import Coordinate
from primerforge.specificity.analyzer import SpecificityAnalyzer
from primerforge.specificity.models import SpecificityResult
from primerforge.specificity.pair_validator import (
    PairSpecificityValidator,
)


class SpecificityEngine:
    """
    Evaluate complete primer pairs for species specificity.
    """

    def __init__(
        self,
        database: str,
        target_species: str,
    ):

        self.database = database

        self.target_species = target_species

        self.analyzer = SpecificityAnalyzer()

        self.pair_validator = PairSpecificityValidator()

    def evaluate_pair(
        self,
        pair: PrimerPair,
        output_dir: Path,
    ) -> SpecificityResult:

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        #
        # Forward primer
        #
        forward_result = self.analyzer.analyse(
            primer=pair.forward,
            database=self.database,
            blast_output=output_dir / "forward.tsv",
            target_species=self.target_species,
        )

        #
        # Reverse primer
        #
        reverse_result = self.analyzer.analyse(
            primer=pair.reverse,
            database=self.database,
            blast_output=output_dir / "reverse.tsv",
            target_species=self.target_species,
        )

        #
        # Pair coordinate
        #
        pair.coordinate = Coordinate(
            start=pair.forward.coordinate.start,
            end=pair.reverse.coordinate.end,
            system=pair.forward.coordinate.system,
        )

        #
        # Combine primer results
        #
        specificity = min(
            forward_result.specificity_score,
            reverse_result.specificity_score,
        )

        passed = (
            forward_result.passed
            and
            reverse_result.passed
        )

        result = SpecificityResult(

            forward_hits=forward_result.forward_hits,

            reverse_hits=reverse_result.forward_hits,

            target_hits=(
                forward_result.target_hits
                +
                reverse_result.target_hits
            ),

            off_target_hits=(
                forward_result.off_target_hits
                +
                reverse_result.off_target_hits
            ),

            specificity_score=specificity,

            passed=passed,

            rejection_reason=(
                None
                if passed
                else "Primer failed specificity"
            ),

        )

        pair.specificity_result = result

        pair.specificity_score = specificity

        pair.passed_specificity = passed

        #
        # Pair-level validation
        #
        validation = self.pair_validator.validate(
            pair,
        )

        if not validation["passed"]:

            pair.passed_specificity = False

            pair.specificity_score = 0.0

            result.passed = False

            result.rejection_reason = (
                "Potential off-target PCR amplification"
            )

        return result
