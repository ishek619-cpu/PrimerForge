"""
Tests for the Coordinate model.
"""

from primerforge.reference.coordinates import Coordinate


def test_coordinate_length():

    coordinate = Coordinate(
        start=100,
        end=120,
    )

    assert coordinate.length == 21


def test_contains():

    coordinate = Coordinate(
        start=10,
        end=20,
    )

    assert coordinate.contains(10)

    assert coordinate.contains(15)

    assert coordinate.contains(20)

    assert not coordinate.contains(9)

    assert not coordinate.contains(21)


def test_overlap():

    a = Coordinate(
        start=10,
        end=20,
    )

    b = Coordinate(
        start=15,
        end=30,
    )

    c = Coordinate(
        start=21,
        end=40,
    )

    assert a.overlaps(b)

    assert not a.overlaps(c)


def test_different_coordinate_systems():

    reference = Coordinate(
        start=10,
        end=20,
        system="reference",
    )

    alignment = Coordinate(
        start=10,
        end=20,
        system="alignment",
    )

    assert not reference.overlaps(
        alignment,
    )


def test_shift():

    coordinate = Coordinate(
        start=100,
        end=120,
    )

    shifted = coordinate.shift(
        50,
    )

    assert shifted.start == 150

    assert shifted.end == 170

    assert shifted.length == 21


def test_to_dict():

    coordinate = Coordinate(
        start=1,
        end=10,
        system="reference",
    )

    data = coordinate.to_dict()

    assert data["start"] == 1

    assert data["end"] == 10

    assert data["length"] == 10

    assert data["system"] == "reference"


def test_repr():

    coordinate = Coordinate(
        start=25,
        end=50,
    )

    text = repr(
        coordinate,
    )

    assert "Coordinate" in text

    assert "25-50" in text
