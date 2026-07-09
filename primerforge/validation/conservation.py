"""
Primer conservation validation.
"""

from pathlib import Path

from Bio import AlignIO

from primerforge.models.primer import Primer


class PrimerConservation:

    """
    Evaluate how well a primer is conserved across an alignment.
    """

    def evaluate(
        self,
        primer: Primer,
        alignment: Path,
        mismatches: int = 1,
    ):

        aln = AlignIO.read(
            alignment,
            "fasta",
        )

        primer_seq = primer.sequence.upper()

        length = len(primer_seq)

        matches = 0

        for record in aln:

            sequence = str(record.seq).upper()

            found = False

            for i in range(
                len(sequence) - length + 1
            ):

                target = sequence[i:i + length]

                diff = sum(
                    a != b
                    for a, b in zip(
                        primer_seq,
                        target,
                    )
                )

                if diff <= mismatches:

                    found = True

                    break

            if found:

                matches += 1

        return {

            "matches": matches,

            "total": len(aln),

            "coverage": round(
                matches / len(aln) * 100,
                2,
            ),

        }
