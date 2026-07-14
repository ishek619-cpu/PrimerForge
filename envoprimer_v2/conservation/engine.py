"""
Population conservation engine.

Calculates primer conservation across the aligned target
population.
"""

from __future__ import annotations

from Bio.Seq import Seq


class ConservationEngine:
    """
    Calculates population conservation for individual primers.
    """

    def __init__(
        self,
        ignore_n: bool = True,
        ignore_gaps: bool = True,
    ):

        self.ignore_n = ignore_n
        self.ignore_gaps = ignore_gaps

    def score_primer(
        self,
        primer,
        alignment: list[str],
    ) -> float:
        """
        Calculate percentage conservation of a primer across an
        aligned target population.

        Parameters
        ----------
        primer
            Primer object.

        alignment
            List of aligned target sequences.

        Returns
        -------
        float
            Conservation percentage.
        """

        if not alignment:
            return 0.0

        #
        # Reverse primer binds to reverse complement.
        #
        sequence = primer.sequence.upper()

        if primer.strand == "-":
            sequence = str(
                Seq(sequence).reverse_complement()
            )

        matches = 0
        total = 0

        for target in alignment:

            target = target.upper()

            if primer.end > len(target):
                continue

            region = target[
                primer.start - 1 : primer.end
            ]

            if len(region) != len(sequence):
                continue

            valid = True

            for a, b in zip(sequence, region):

                if self.ignore_gaps and b == "-":
                    valid = False
                    break

                if self.ignore_n and b == "N":
                    valid = False
                    break

                if a != b:
                    valid = False
                    break

            total += 1

            if valid:
                matches += 1

        if total == 0:
            return 0.0

        return round(
            matches * 100 / total,
            2,
        )
