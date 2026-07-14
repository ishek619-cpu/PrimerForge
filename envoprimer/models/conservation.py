"""
Conservation profile model.
"""

from dataclasses import dataclass


@dataclass(slots=True)
class ConservationProfile:
    """
    Conservation scores for an alignment.
    """

    alignment_length: int

    number_of_sequences: int

    scores: list[float]

    @property
    def average(self) -> float:

        if not self.scores:
            return 0.0

        return sum(self.scores) / len(self.scores)

    def score(self, position: int) -> float:

        return self.scores[position]
