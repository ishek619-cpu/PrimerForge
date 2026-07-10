"""
Primer conservation analysis.
"""

from pathlib import Path

from Bio import AlignIO


class PrimerConservation:
    """
    Evaluate primer conservation from a multiple sequence alignment.
    """

    def __init__(self):

        self._alignment = None

        self._alignment_path = None

    def _load_alignment(
        self,
        alignment: Path,
    ):

        if (
            self._alignment is None
            or self._alignment_path != alignment
        ):

            self._alignment = AlignIO.read(
                alignment,
                "fasta",
            )

            self._alignment_path = alignment

        return self._alignment

    def evaluate(
        self,
        primer,
        alignment: Path,
    ) -> dict:

        aln = self._load_alignment(
            alignment,
        )

        length = len(primer.sequence)

        best_score = -1.0
        best_start = -1

        conservation = []

        for start in range(
            aln.get_alignment_length() - length + 1
        ):

            matches = 0
            total = 0

            for i in range(length):

                column = aln[:, start + i]

                bases = [
                    b
                    for b in column
                    if b != "-"
                ]

                if not bases:
                    continue

                total += 1

                most_common = max(
                    set(bases),
                    key=bases.count,
                )

                freq = (
                    bases.count(most_common)
                    / len(bases)
                )

                matches += freq

            score = (
                matches / total * 100
                if total
                else 0.0
            )

            conservation.append(
                round(score, 2)
            )

            if score > best_score:

                best_score = score

                best_start = start

        return {

            "best_start": best_start,

            "window_length": length,

            "best_score": round(
                best_score,
                2,
            ),

            "profile": conservation,

        }
