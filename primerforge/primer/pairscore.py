"""
Primer pair scoring engine.
"""

from primerforge.models.pair import PrimerPair


class PrimerPairScorer:
    """
    Score and rank primer pairs.
    """

    def score(self, pair: PrimerPair) -> PrimerPair:

        score = 100.0

        # Tm difference
        dtm = abs(pair.forward.tm - pair.reverse.tm)
        score -= dtm * 5

        # GC difference
        dgc = abs(pair.forward.gc - pair.reverse.gc)
        score -= dgc * 0.5

        # Product size
        if not (80 <= pair.product_size <= 250):
            score -= 20

        # GC clamps
        if pair.forward.gc_clamp:
            score += 2

        if pair.reverse.gc_clamp:
            score += 2

        # Hairpins
        score -= pair.forward.hairpin_score
        score -= pair.reverse.hairpin_score

        # Self dimers
        score -= pair.forward.self_dimer_score
        score -= pair.reverse.self_dimer_score

        pair.score = round(score, 2)

        return pair

    def rank(self, pairs):

        return sorted(
            [self.score(p) for p in pairs],
            key=lambda p: p.score,
            reverse=True,
        )
