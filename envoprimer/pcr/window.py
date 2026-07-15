"""
Window extraction for seed-and-extend PCR alignment.
"""

from __future__ import annotations


class WindowExtractor:

    def __init__(
        self,
        flank: int = 25,
    ):

        self.flank = flank

    ###########################################################

    def extract(
        self,
        sequence: str,
        seed_position: int,
        seed_length: int,
    ):

        start = max(
            0,
            seed_position - self.flank,
        )

        end = min(
            len(sequence),
            seed_position + seed_length + self.flank,
        )

        return sequence[start:end], start
