"""
Tests for SNP-aware primer optimization.
"""

from envoprimer.models.primer import Primer
from envoprimer.models.snp import SNP
from envoprimer.primer.snp_optimizer import SNPOptimizer


def test_snp_optimizer():

    primer = Primer(
        sequence="ATGCGTACGTAGCTAGCTAG",
        start=100,
        end=119,
        strand="+",
        length=20,
    )

    snps = [
        SNP(
            position=118,
            reference="A",
            alternatives=["G"],
            counts={
                "A": 20,
                "G": 5,
            },
        )
    ]

    score = SNPOptimizer().score(
        primer,
        snps,
    )

    assert score < 100
