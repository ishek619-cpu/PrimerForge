"""
Advanced BLAST specificity scoring.
"""


class SpecificityScorer:
    """
    Score primer specificity using BLAST hits.
    """

    def score(
        self,
        hits,
    ) -> float:

        if not hits:
            return 100.0

        score = 100.0

        perfect_hits = 0

        for hit in hits:

            identity = hit.get(
                "identity",
                0.0,
            )

            mismatches = hit.get(
                "mismatches",
                0,
            )

            length = hit.get(
                "length",
                0,
            )

            bitscore = hit.get(
                "bitscore",
                0.0,
            )

            evalue = hit.get(
                "evalue",
                1.0,
            )

            # Ignore weak alignments
            if evalue > 1e-3:
                continue

            # Perfect full-length match
            if (
                identity == 100.0
                and mismatches == 0
            ):

                perfect_hits += 1

                if perfect_hits > 1:
                    score -= 20.0

                continue

            # Near-perfect off-target
            if identity >= 99.0:

                score -= 10.0

            elif identity >= 95.0:

                score -= 6.0

            elif identity >= 90.0:

                score -= 3.0

            # Long alignments are more concerning
            if length >= 18:

                score -= 2.0

            elif length >= 15:

                score -= 1.0

            # Strong BLAST hit
            if bitscore >= 35:

                score -= 1.0

        return round(
            max(score, 0.0),
            2,
        )
