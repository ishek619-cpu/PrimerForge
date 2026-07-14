"""
Tests for multiplex primer compatibility.
"""

from envoprimer.models.primer import Primer
from envoprimer.models.pair import PrimerPair

from envoprimer.primer.multiplex import (
    MultiplexAnalyzer,
)


class DummyThermo:

    def heterodimer(
        self,
        primer1,
        primer2,
    ):

        return 5.0


def test_multiplex():

    p1 = PrimerPair(

        forward=Primer(
            sequence="AAAAAAAAAAAAAAAAAAAA",
            start=1,
            end=20,
            strand="+",
            length=20,
        ),

        reverse=Primer(
            sequence="CCCCCCCCCCCCCCCCCCCC",
            start=101,
            end=120,
            strand="-",
            length=20,
        ),

        product_size=120,

    )

    p2 = PrimerPair(

        forward=Primer(
            sequence="GGGGGGGGGGGGGGGGGGGG",
            start=201,
            end=220,
            strand="+",
            length=20,
        ),

        reverse=Primer(
            sequence="TTTTTTTTTTTTTTTTTTTT",
            start=301,
            end=320,
            strand="-",
            length=20,
        ),

        product_size=120,

    )

    analyzer = MultiplexAnalyzer()

    analyzer.thermo = DummyThermo()

    result = analyzer.evaluate(
        [
            p1,
            p2,
        ]
    )

    assert result["compatible_pairs"] == 1

    assert result["conflicting_pairs"] == 0

    assert len(result["compatible"]) == 1
