"""
Tests for primer pair specificity.
"""

from primerforge.models.primer import Primer
from primerforge.models.pair import PrimerPair

from primerforge.specificity.pair_analyzer import (
    PrimerPairSpecificity,
)


class DummyScorer:

    def score(self, hits):

        return hits[0]["score"]


class DummyValidator:

    def __init__(self):

        self.scorer = DummyScorer()

    def blast_hits(
        self,
        sequence,
        database,
    ):

        if sequence.startswith("AAA"):

            return [

                {
                    "score": 100.0,
                    "subject": "chr1",
                    "query_start": 100,
                    "query_end": 120,
                    "identity": 100.0,
                }

            ]

        return [

            {
                "score": 90.0,
                "subject": "chr1",
                "query_start": 250,
                "query_end": 270,
                "identity": 90.0,
            }

        ]


def test_pair_specificity():

    forward = Primer(
        sequence="AAAAAAAAAAAAAAAAAAAA",
        start=1,
        end=20,
        strand="+",
        length=20,
    )

    reverse = Primer(
        sequence="CCCCCCCCCCCCCCCCCCCC",
        start=100,
        end=119,
        strand="-",
        length=20,
    )

    pair = PrimerPair(
        forward=forward,
        reverse=reverse,
        product_size=120,
        score=80.0,
    )

    analyzer = PrimerPairSpecificity()

    analyzer.validator = DummyValidator()

    result = analyzer.evaluate(
        pair,
        "dummy",
    )

    assert result["forward_score"] == 100.0

    assert result["reverse_score"] == 90.0

    assert len(result["products"]) == 1

    assert result["pair_score"] == 95.0

    assert result["passed"] is True

    pair = analyzer.annotate(
        pair,
        "dummy",
    )

    assert pair.score == 76.0
