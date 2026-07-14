"""
Off-target amplification risk assessment.
"""

from envoprimer.models.pcr_product import PCRProduct
from envoprimer.specificity.models import BlastHit


class RiskAssessor:
    """
    Estimate the probability that an off-target BLAST hit
    could generate a PCR product.

    Final score is between 0 and 100.

    Lower score = safer primer.
    """

    def assess(
        self,
        hit: BlastHit,
        product: PCRProduct | None = None,
    ) -> dict:

        score = 0.0

        #
        # Sequence identity
        #
        if hit.identity >= 99:

            score += 35

        elif hit.identity >= 97:

            score += 30

        elif hit.identity >= 95:

            score += 25

        elif hit.identity >= 90:

            score += 15

        else:

            score += 5

        #
        # Query coverage
        #
        if hit.coverage >= 100:

            score += 20

        elif hit.coverage >= 95:

            score += 15

        elif hit.coverage >= 90:

            score += 10

        #
        # 3' mismatch score
        #
        if hit.three_prime_score >= 95:

            score += 20

        elif hit.three_prime_score >= 80:

            score += 15

        elif hit.three_prime_score >= 60:

            score += 8

        else:

            score += 0

        #
        # Number of 3' mismatches
        #
        if hit.three_prime_mismatches == 0:

            score += 10

        elif hit.three_prime_mismatches == 1:

            score += 5

        #
        # Alignment quality
        #
        if hit.bitscore >= 100:

            score += 10

        elif hit.bitscore >= 70:

            score += 6

        elif hit.bitscore >= 50:

            score += 3

        #
        # E-value
        #
        if hit.evalue <= 1e-20:

            score += 5

        elif hit.evalue <= 1e-10:

            score += 3

        elif hit.evalue <= 1e-5:

            score += 1

        #
        # Predicted PCR product
        #
        if product is not None:

            if product.passed:

                score += 20

            else:

                score += 5

        #
        # Clamp
        #
        score = min(
            round(score, 2),
            100.0,
        )

        #
        # Risk category
        #
        if score >= 80:

            level = "CRITICAL"

        elif score >= 60:

            level = "HIGH"

        elif score >= 40:

            level = "MODERATE"

        else:

            level = "LOW"

        return {

            "score": score,

            "level": level,

            "should_reject": level in {
                "HIGH",
                "CRITICAL",
            },

        }
