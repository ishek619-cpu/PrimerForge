"""
Population-wide primer evaluation.
"""

from pathlib import Path

from primerforge.analysis.conservation import (
    ConservedRegionFinder,
)


class PopulationAnalyzer:
    """
    Evaluate primer robustness across all available
    sequences of the target species.
    """

    def __init__(self):

        self.conservation = ConservedRegionFinder()

    def analyse(
        self,
        primer: str,
        alignment: Path,
        start: int,
    ) -> dict:

        conservation = self.conservation.primer_conservation(
            primer=primer,
            alignment=alignment,
            start=start,
        )

        return {

            "population_conservation": conservation,

            "population_coverage": conservation,

            "estimated_failure_rate": round(
                100.0 - conservation,
                2,
            ),

        }
