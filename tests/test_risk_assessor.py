"""
Tests for off-target risk assessment.
"""

from envoprimer.models.pcr_product import PCRProduct
from envoprimer.risk.assessor import RiskAssessor
from envoprimer.specificity.models import BlastHit


def test_high_risk():

    hit = BlastHit(
        accession="NC_000001",
        species="Species A",
        identity=99.8,
        coverage=100.0,
        alignment_length=25,
        mismatches=0,
        gap_opens=0,
        qstart=1,
        qend=25,
        sstart=100,
        send=124,
        query_sequence="ATGCGTACGTAGCTAGCTAGCTAGC",
        subject_sequence="ATGCGTACGTAGCTAGCTAGCTAGC",
        strand="plus",
        bitscore=120.0,
        evalue=1e-40,
        three_prime_mismatches=0,
        three_prime_score=100.0,
    )

    product = PCRProduct(
        chromosome="NC_000001",
        forward_start=100,
        reverse_end=250,
        product_size=151,
        forward_strand="+",
        reverse_strand="-",
        passed=True,
    )

    result = RiskAssessor().assess(
        hit,
        product,
    )

    assert result["level"] == "CRITICAL"

    assert result["should_reject"] is True

    assert result["score"] >= 80


def test_low_risk():

    hit = BlastHit(
        accession="NC_000002",
        species="Species B",
        identity=84.0,
        coverage=80.0,
        alignment_length=18,
        mismatches=6,
        gap_opens=1,
        qstart=1,
        qend=18,
        sstart=500,
        send=517,
        query_sequence="ATGCGTACGTAGCTAGCT",
        subject_sequence="ATACGTTCGTGGCTTGCT",
        strand="plus",
        bitscore=32.0,
        evalue=0.1,
        three_prime_mismatches=3,
        three_prime_score=30.0,
    )

    result = RiskAssessor().assess(hit)

    assert result["level"] == "LOW"

    assert result["should_reject"] is False

    assert result["score"] < 40
