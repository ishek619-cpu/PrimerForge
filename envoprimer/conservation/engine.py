"""
Population conservation engine.

Calculates how well each primer is conserved across all target
sequences.
"""

from __future__ import annotations

from Bio.Seq import Seq


class ConservationEngine:

    def __init__(

        self,

        ignore_n: bool = True,

        ignore_gap: bool = True,

    ):

        self.ignore_n = ignore_n

        self.ignore_gap = ignore_gap

    def score(

        self,

        pair,

        alignment: list[str],

    ):

        pair.forward.population_score = self.score_primer(

            pair.forward,

            alignment,

        )

        pair.reverse.population_score = self.score_primer(

            pair.reverse,

            alignment,

        )

        pair.population_score = min(

            pair.forward.population_score,

            pair.reverse.population_score,

        )

        return pair.population_score

    def score_primer(

        self,

        primer,

        alignment,

    ):

        if not alignment:

            return 0.0

        sequence = primer.sequence.upper()

        if primer.strand == "-":

            sequence = str(

                Seq(sequence).reverse_complement()

            )

        total = 0

        perfect = 0

        for target in alignment:

            target = target.upper()

            if primer.end > len(target):

                continue

            region = target[

                primer.start - 1:

                primer.end

            ]

            if len(region) != len(sequence):

                continue

            mismatch = False

            for a, b in zip(

                sequence,

                region,

            ):

                if self.ignore_gap and b == "-":

                    mismatch = True

                    break

                if self.ignore_n and b == "N":

                    mismatch = True

                    break

                if a != b:

                    mismatch = True

                    break

            total += 1

            if not mismatch:

                perfect += 1

        if total == 0:

            return 0.0

        return round(

            perfect * 100.0 / total,

            2,

        )
