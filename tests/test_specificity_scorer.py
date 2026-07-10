"""
Tests for BLAST specificity scoring.
"""

from primerforge.specificity.models import BlastHit
from primerforge.specificity.scorer import SpecificityScorer


def test_specificity_score():

    scorer = SpecificityScorer()

    hits = [

        BlastHit(

            accession="NC_013663.1",

            species="Oreochromis niloticus",

            identity=100.0,

            coverage=100.0,

            alignment_length=20,

            mismatches=0,

            gap_opens=0,

            qstart=1,

            qend=20,

            sstart=100,

            send=119,

            strand="plus",

            bitscore=40.0,

            evalue=1e-10,

        )

    ]

    score = scorer.score(hits)

    assert score == 100.0
