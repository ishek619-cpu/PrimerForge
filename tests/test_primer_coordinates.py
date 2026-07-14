"""
Tests for Primer coordinate integration.
"""

from envoprimer.models.primer import Primer
from envoprimer.reference.coordinates import Coordinate


def test_default_coordinate():

    primer = Primer(
        sequence="ATGCATGC",
        start=101,
        end=108,
        strand="+",
        length=8,
    )

    assert primer.coordinate.start == 101

    assert primer.coordinate.end == 108

    assert primer.coordinate.system == "reference"


def test_custom_coordinate():

    coordinate = Coordinate(
        start=500,
        end=507,
        system="alignment",
    )

    primer = Primer(
        sequence="ATGCATGC",
        start=1,
        end=8,
        strand="+",
        length=8,
        coordinate=coordinate,
    )

    assert primer.coordinate.system == "alignment"

    assert primer.coordinate.start == 500


def test_coordinate_length():

    primer = Primer(
        sequence="ATGCATGC",
        start=25,
        end=32,
        strand="+",
        length=8,
    )

    assert primer.coordinate_length == 8


def test_to_dict():

    primer = Primer(
        sequence="ATGCATGC",
        start=10,
        end=17,
        strand="+",
        length=8,
        tm=60.2,
        gc=50.0,
    )

    data = primer.to_dict()

    assert data["sequence"] == "ATGCATGC"

    assert data["coordinate"]["start"] == 10

    assert data["coordinate"]["end"] == 17

    assert data["coordinate"]["system"] == "reference"


def test_repr():

    primer = Primer(
        sequence="ATGC",
        start=1,
        end=4,
        strand="+",
        length=4,
    )

    text = repr(primer)

    assert "Primer(" in text

    assert "Coordinate" in text
