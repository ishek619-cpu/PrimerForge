"""
Alignment reconstruction utilities.

These utilities reconstruct the full primer alignment
using BLAST coordinates so downstream analyses operate
on primer coordinates instead of BLAST alignment
coordinates.
"""

from dataclasses import dataclass

from primerforge.specificity.models import BlastHit


@dataclass(slots=True)
class AlignmentPosition:
    """
    One position within the primer.
    """

    primer_position: int

    primer_base: str

    subject_base: str

    aligned: bool

    mismatch: bool


class AlignmentReconstructor:
    """
    Convert BLAST alignments into a full primer alignment.
    """

    def reconstruct(
        self,
        hit: BlastHit,
        primer_length: int,
    ) -> list[AlignmentPosition]:

        positions = []

        #
        # Initially every base is unaligned.
        #
        for i in range(primer_length):

            positions.append(

                AlignmentPosition(

                    primer_position=i + 1,

                    primer_base="N",

                    subject_base="-",

                    aligned=False,

                    mismatch=False,

                )

            )

        query = hit.query_sequence.upper()

        subject = hit.subject_sequence.upper()

        primer_index = hit.qstart - 1

        for q, s in zip(query, subject):

            #
            # Skip gaps in the primer.
            #
            if q == "-":
                continue

            mismatch = (
                q != s
                and s != "-"
            )

            positions[primer_index] = AlignmentPosition(

                primer_position=primer_index + 1,

                primer_base=q,

                subject_base=s,

                aligned=True,

                mismatch=mismatch,

            )

            primer_index += 1

        return positions
