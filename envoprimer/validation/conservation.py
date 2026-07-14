"""
Primer conservation analysis.
"""

from pathlib import Path

from Bio import AlignIO


class PrimerConservation:
    """
    Evaluate conservation only at the designed primer position.
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

        #
        # Primer coordinate (1-based -> 0-based)
        #
        start = primer.coordinate.start - 1

        length = len(
            primer.sequence,
        )

        alignment_length = aln.get_alignment_length()

        if start < 0 or start + length > alignment_length:

            return {
                "best_start": start,
                "window_length": length,
                "best_score": 0.0,
                "profile": [],
            }

        matches = 0.0

        total = 0

        profile = []

        for i in range(length):

            column = aln[:, start + i]

            bases = [
                b
                for b in column
                if b != "-"
            ]

            if not bases:

                profile.append(0.0)

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

            profile.append(
                round(
                    freq * 100,
                    2,
                )
            )

            matches += freq

        score = (
            matches
            / total
            * 100
            if total
            else 0.0
        )

        return {

            "best_start": start,

            "window_length": length,

            "best_score": round(
                score,
                2,
            ),

            "profile": profile,

        }
