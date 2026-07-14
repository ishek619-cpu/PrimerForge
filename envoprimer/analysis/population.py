"""
Population-wide primer evaluation.
"""

from pathlib import Path

from envoprimer.analysis.conservation import ConservedRegionFinder
from envoprimer.reference.coordinates import Coordinate


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
        start: int | None = None,
        coordinate: Coordinate | None = None,
    ) -> dict:

        if coordinate is not None:
            start = coordinate.start

        if start is None:
            raise ValueError(
                "Either start or coordinate must be provided."
            )

        conservation = self.conservation.primer_conservation(
            primer=primer,
            alignment=alignment,
            start=start,
        )

        print("\n========== POPULATION ==========")
        print("Primer:", primer)
        print("Start :", start)
        print("Conservation:", conservation)
        print("Stats:", self.conservation.last_statistics)
        print("================================\n")

        return {

            "population_conservation": conservation,

            "population_coverage": conservation,

            "estimated_failure_rate": round(
                100.0 - conservation,
                2,
            ),

        }
