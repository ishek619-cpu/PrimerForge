"""
Forward-strand seed search for PCR.

Convention:
- Primers are always supplied exactly as written (5'→3').
- Forward primer searches use the primer sequence.
- Reverse primer searches use the reverse complement.
"""

from __future__ import annotations

from collections import defaultdict


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


class SeedSearcher:
    """
    K-mer seed index built on the forward genomic strand.

    All coordinates returned are in the original genome coordinate system.
    """

    def __init__(self, seed_length: int = 18):

        self.seed_length = seed_length
        self.index = defaultdict(list)

    def build(self, sequences: list[str]) -> None:
        """Build a forward-strand seed index."""

        self.index.clear()

        for seq_id, sequence in enumerate(sequences):

            sequence = sequence.upper()

            for i in range(len(sequence) - self.seed_length + 1):

                seed = sequence[i:i + self.seed_length]

                self.index[seed].append((seq_id, i))

    def search(
        self,
        primer: str,
        *,
        reverse: bool = False,
    ):
        """
        Search for primer seed.

        Parameters
        ----------
        primer
            Primer sequence (5'→3').

        reverse
            False → search primer directly.
            True  → search reverse-complement of primer.
        """

        primer = primer.upper()

        if reverse:
            primer = reverse_complement(primer)

        seed = primer[-self.seed_length:]

        return list(self.index.get(seed, []))
