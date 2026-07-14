"""
Tests for nearest-neighbor thermodynamic calculations.
"""

from envoprimer.thermodynamics.nearest_neighbor import (
    NearestNeighborCalculator,
)


def test_calculate():

    calculator = NearestNeighborCalculator()

    result = calculator.calculate(
        "ATGCGTACGA",
    )

    assert result.delta_h < 0

    assert result.delta_s < 0

    assert isinstance(
        result.delta_g,
        float,
    )

    assert isinstance(
        result.tm,
        float,
    )


def test_short_sequence():

    calculator = NearestNeighborCalculator()

    result = calculator.calculate(
        "AT",
    )

    assert result.delta_h < 0

    assert result.delta_s < 0


def test_unknown_bases():

    calculator = NearestNeighborCalculator()

    result = calculator.calculate(
        "ATNNGC",
    )

    assert isinstance(
        result.tm,
        float,
    )

    assert isinstance(
        result.delta_g,
        float,
    )
