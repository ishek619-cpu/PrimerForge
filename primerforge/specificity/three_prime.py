"""
Exact 3' mismatch analysis using reconstructed primer alignments.
"""

from primerforge.specificity.alignment import (
    AlignmentReconstructor,
)
from primerforge.specificity.models import BlastHit
from primerforge.specificity.mismatch_matrix import (
    mismatch_penalty,
)


class ThreePrimeAnalyzer:
    """
    Analyse the last N primer bases.
    """

    def __init__(
        self,
        window: int = 5,
    ):

        self.window = window

        self.reconstructor = AlignmentReconstructor()

    def mismatches(
        self,
        hit: BlastHit,
    ) -> list[dict]:

        primer_length = max(
            len(hit.query_sequence),
            hit.qend,
        )

        alignment = self.reconstructor.reconstruct(
            hit,
            primer_length,
        )

        region = alignment[-self.window:]

        mismatches = []

        for position in region:

            #
            # Unaligned bases count as unknown
            #
            if not position.aligned:

                mismatches.append(
                    {
                        "position": position.primer_position,
                        "query": position.primer_base,
                        "subject": "-",
                        "aligned": False,
                    }
                )

                continue

            if position.mismatch:

                mismatches.append(
                    {
                        "position": position.primer_position,
                        "query": position.primer_base,
                        "subject": position.subject_base,
                        "aligned": True,
                    }
                )

        return mismatches

    def count(
        self,
        hit: BlastHit,
    ) -> int:

        return len(
            self.mismatches(
                hit,
            )
        )

    def score(
        self,
        hit: BlastHit,
    ) -> float:

        mismatches = self.mismatches(
            hit,
        )

        if not mismatches:
            return 100.0

        primer_length = max(
            len(hit.query_sequence),
            hit.qend,
        )

        score = 100.0

        for mismatch in mismatches:

            #
            # Unaligned terminal bases
            #
            if not mismatch["aligned"]:

                score -= 25

                continue

            chemistry = mismatch_penalty(
                mismatch["query"],
                mismatch["subject"],
            )

            distance = (
                primer_length
                - mismatch["position"]
            )

            #
            # Position weighting
            #
            if distance == 0:

                weight = 40

            elif distance == 1:

                weight = 30

            elif distance == 2:

                weight = 20

            elif distance == 3:

                weight = 10

            else:

                weight = 5

            score -= chemistry * weight

        return max(
            round(score, 2),
            0.0,
        )
