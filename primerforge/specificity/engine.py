"""
Species-specific primer evaluation engine.
"""

from pathlib import Path

from primerforge.models.pair import PrimerPair

from primerforge.specificity.analyzer import SpecificityAnalyzer
from primerforge.specificity.models import SpecificityResult


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

    def evaluate_pair(
        self,
        pair: PrimerPair,
        output_dir: Path,
    ) -> SpecificityResult:

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        forward_result = self.analyzer.analyse(
            primer=pair.forward,
            database=self.database,
            blast_output=output_dir / "forward.tsv",
            target_species=self.target_species,
        )

        reverse_result = self.analyzer.analyse(
            primer=pair.reverse,
            database=self.database,
            blast_output=output_dir / "reverse.tsv",
            target_species=self.target_species,
        )

        specificity = min(
            forward_result.specificity_score,
            reverse_result.specificity_score,
        )

        passed = (
            forward_result.passed
            and reverse_result.passed
        )

        result = SpecificityResult(
            forward_hits=forward_result.forward_hits,
            reverse_hits=reverse_result.forward_hits,
            target_hits=(
                forward_result.target_hits
                + reverse_result.target_hits
            ),
            off_target_hits=(
                forward_result.off_target_hits
                + reverse_result.off_target_hits
            ),
            specificity_score=specificity,
            passed=passed,
            rejection_reason=(
                None
                if passed
                else "Primer pair failed specificity"
            ),
        )

        pair.specificity_score = specificity
        pair.passed_specificity = passed
        pair.specificity_result = result

        return result
