"""
Tests for primer pair specificity.
"""

from primerforge.models.primer import Primer
from primerforge.models.pair import PrimerPair

from primerforge.specificity.pair_analyzer import (
    PrimerPairSpecificity,
)


class DummyValidator:

    def validate(
        self,
        sequence,
        database,
    ):

        if sequence.startswith("AAA"):
            return 100.0

        return 90.0


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

    assert result["pair_score"] == 95.0

    assert result["passed"] is True

    pair = analyzer.annotate(
        pair,
        "dummy",
    )

    assert pair.score == 76.0

    result = analyzer.evaluate(
        pair,
        "dummy",
    )

    assert result["passed"] is True
