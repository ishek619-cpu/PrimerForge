"""
Primer conservation analysis.
"""

from pathlib import Path

from Bio import AlignIO

from primerforge.models.primer import Primer


class PrimerConservation:
    """
    Analyse primer conservation across an alignment.
    """

    def evaluate(
        self,
        primer: Primer,
        alignment: Path,
    ):

        aln = AlignIO.read(
            alignment,
            "fasta",
        )

        sequence = primer.sequence.upper()

        length = len(sequence)

        exact = 0
        one = 0
        two = 0

        mismatches = []

        for record in aln:

            target = str(
                record.seq[
                    primer.start:
                    primer.start + length
                ]
            ).upper()

            diff = sum(
                a != b
                for a, b in zip(
                    sequence,
                    target,
                )
            )

            mismatches.append(diff)

            if diff == 0:

                exact += 1

            elif diff == 1:

                one += 1

            elif diff == 2:

                two += 1

        coverage = (
            (exact + one)
            / len(aln)
            * 100
        )

        return {

            "exact": exact,

            "one_mismatch": one,

            "two_mismatch": two,

            "coverage": round(
                coverage,
                2,
            ),

            "mean_mismatches": round(
                sum(mismatches)
                / len(mismatches),
                2,
            ),

        }
