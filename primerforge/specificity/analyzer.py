"""
Primer specificity analysis.
"""

from primerforge.specificity.blast import BlastRunner


class SpecificityAnalyzer:
    """
    Analyse BLAST hits for primer specificity.
    """

    def __init__(self):

        self.blast = BlastRunner()

    def analyse(
        self,
        primer,
        database,
    ):

        hits = self.blast.search(
            primer.sequence,
            database,
        )

        perfect = 0
        offtargets = 0
        best_identity = 0.0
        best_coverage = 0.0

        for hit in hits:

            identity = hit.get(
                "identity",
                0.0,
            )

            coverage = hit.get(
                "coverage",
                0.0,
            )

            if identity > best_identity:

                best_identity = identity

            if coverage > best_coverage:

                best_coverage = coverage

            if identity == 100.0:

                perfect += 1

            elif identity >= 90.0:

                offtargets += 1

        score = 100.0

        score -= offtargets * 5

        if perfect > 1:

            score -= (
                perfect - 1
            ) * 10

        score = max(
            score,
            0,
        )

        return {

            "hits": len(hits),

            "perfect_hits": perfect,

            "offtargets": offtargets,

            "best_identity": round(
                best_identity,
                2,
            ),

            "best_coverage": round(
                best_coverage,
                2,
            ),

            "score": round(
                score,
                2,
            ),

        }
