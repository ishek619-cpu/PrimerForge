"""
Coordinate search window utilities.
"""

from __future__ import annotations


class CoordinateWindow:
    """
    Generate candidate alignment positions around the
    expected primer coordinate.

    This replaces exhaustive scanning of every alignment
    position.
    """

    def __init__(
        self,
        tolerance: int = 3,
    ):
        self.tolerance = tolerance

    def starts(
        self,
        expected: int,
        sequence_length: int,
    ) -> range:
        """
        Return candidate alignment starts.

        Coordinates entering this function are 1-based.
        PrimerMatcher expects 0-based indices.
        """

        expected -= 1

        start = max(
            0,
            expected - self.tolerance,
        )

        end = min(
            sequence_length - 1,
            expected + self.tolerance,
        )

        return range(
            start,
            end + 1,
        )
