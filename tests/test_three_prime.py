"""
Tests for exact 3' mismatch analysis.
"""

from primerforge.specificity.models import BlastHit
from primerforge.specificity.three_prime import ThreePrimeAnalyzer


def test_three_prime_analysis():

    hit = BlastHit(

        accession="ABC123",

        species="Oreochromis niloticus",

        identity=95.0,

        coverage=100.0,

        alignment_length=20,

        mismatches=2,

        gap_opens=0,

        qstart=1,

        qend=20,

        sstart=100,

        send=119,

        query_sequence="ATCGATCGATCGATCGATCG",

        subject_sequence="ATCGATCGATCGATCGATGG",

        strand="plus",

        bitscore=40.0,

        evalue=1e-10,

    )

    analyzer = ThreePrimeAnalyzer()

    mismatches = analyzer.mismatches(hit)

    assert len(mismatches) == 1

    assert mismatches[0]["position"] == 19

    assert mismatches[0]["query"] == "C"

    assert mismatches[0]["subject"] == "G"

    assert analyzer.count(hit) == 1

    score = analyzer.score(hit)

    assert score < 100.0

    assert score > 0.0
