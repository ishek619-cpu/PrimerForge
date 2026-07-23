"""
Extract local sequence windows around candidate primer binding sites.
"""

from __future__ import annotations

from typing import Tuple


class WindowExtractor:
    """
    Extract a fixed-size sequence window around a candidate binding site.
    """

    def __init__(self, flank: int = 25):
        self.flank = flank

    def extract(
        self,
        sequence: str,
        seed_start: int,
        seed_length: int,
    ) -> Tuple[str, int]:
        """
        Extract a sequence window centered on the seed.

        Parameters
        ----------
        sequence
            Full DNA sequence.

        seed_start
            Position of seed match.

        seed_length
            Length of matching seed.

        Returns
        -------
        window
            Extracted sequence.

        window_start
            Position of window relative to original sequence.
        """

        left = max(0, seed_start - self.flank)
        right = min(len(sequence), seed_start + seed_length + self.flank)

        return sequence[left:right], left
