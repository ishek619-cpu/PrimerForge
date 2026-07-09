"""
BLAST specificity scoring.
"""


class SpecificityScorer:
    """
    Score primer specificity.
    """

    def score(self, hits):

        if len(hits) == 0:
            return 100.0

        perfect = 0

        for hit in hits:

            if (
                hit["identity"] >= 100.0
                and hit["mismatches"] == 0
            ):
                perfect += 1

        score = max(
            0.0,
            100.0 - (perfect - 1) * 20.0,
        )

        return round(score, 2)
