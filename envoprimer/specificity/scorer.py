"""
Advanced BLAST specificity scoring.
"""

from envoprimer.specificity.models import BlastHit


class SpecificityScorer:
    """
    Score primer specificity using BLAST hits.
    """

    def score(
        self,
        hits: list[BlastHit],
    ) -> float:

        if not hits:
            return 100.0

        score = 100.0

        perfect_hits = 0

        for hit in hits:

            # Ignore weak alignments
            if hit.evalue > 1e-3:
                continue

            # Perfect full-length match
            if (
                hit.identity == 100.0
                and hit.mismatches == 0
            ):

                perfect_hits += 1

                if perfect_hits > 1:
                    score -= 20.0

                continue

            # Near-perfect off-target
            if hit.identity >= 99.0:

                score -= 10.0

            elif hit.identity >= 95.0:

                score -= 6.0

            elif hit.identity >= 90.0:

                score -= 3.0

            # Long alignments are more concerning
            if hit.alignment_length >= 18:

                score -= 2.0

            elif hit.alignment_length >= 15:

                score -= 1.0

            # Strong BLAST hit
            if hit.bitscore >= 35:

                score -= 1.0

            # Penalize 3′ mismatches less heavily
            if hit.three_prime_mismatches > 0:

                score += min(
                    hit.three_prime_mismatches,
                    3,
                )

        return round(
            max(min(score, 100.0), 0.0),
            2,
        )
