"""
PrimerForge Primer Candidate Generator
"""

from primerforge.models.primer import Primer
from primerforge.primer.score import PrimerScorer


class PrimerDesigner:
    """
    Generate candidate primers from DNA sequences.
    """

    def __init__(self):

        self.scorer = PrimerScorer()

    def generate_candidates(
        self,
        sequence: str,
        min_length: int = 20,
        max_length: int = 24,
    ) -> list[Primer]:

        candidates = []

        sequence = sequence.upper()

        for length in range(
            min_length,
            max_length + 1,
        ):

            for start in range(
                len(sequence) - length + 1,
            ):

                fragment = sequence[
                    start:start + length
                ]

                primer = Primer(
                    sequence=fragment,
                    start=start + 1,
                    end=start + length,
                    strand="+",
                    length=length,
                )

                primer = self.scorer.evaluate(
                    primer
                )

                candidates.append(
                    primer
                )

        return candidates
