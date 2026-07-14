"""
Tests for PrimerPairScorer.
"""

from envoprimer.models.primer import Primer
from envoprimer.models.pair import PrimerPair
from envoprimer.primer.pairscore import PrimerPairScorer


def test_pairscore():

    forward = Primer(
        sequence="AAAAAAAAAAAAAAAAAAAA",
        start=1,
        end=20,
        strand="+",
        length=20,
    )

    reverse = Primer(
        sequence="TTTTTTTTTTTTTTTTTTTT",
        start=120,
        end=139,
        strand="-",
        length=20,
    )

    pair = PrimerPair(
        forward=forward,
        reverse=reverse,
        product_size=150,
        score=80.0,
    )

    validation = {
        "thermo": 95.0,
        "conservation": 98.0,
        "snp": 100.0,
    }

    specificity = {
        "specificity": 99.0,
    }

    multiplex = {
        "compatibility": 97.0,
    }

    scorer = PrimerPairScorer()

    score = scorer.score(
        pair,
        validation,
        specificity,
        multiplex,
    )

    assert score > 80

    ranked = scorer.rank(
        [pair]
    )

    assert ranked[0] == pair
