"""
SNP-aware primer optimization.
"""

from envoprimer.models.primer import Primer
from envoprimer.models.snp import SNP


class SNPOptimizer:
    """
    Penalize primers that contain SNPs,
    especially near the 3' end.
    """

    def score(
        self,
        primer: Primer,
        snps: list[SNP],
    ) -> float:

        score = 100.0

        for snp in snps:

            if not (
                primer.start <= snp.position <= primer.end
            ):
                continue

            distance = primer.end - snp.position

            if distance <= 2:
                score -= 50

            elif distance <= 5:
                score -= 25

            else:
                score -= 5

        return max(
            score,
            0.0,
        )
