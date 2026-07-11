"""
Universal coordinate mapper.
"""

from __future__ import annotations

from Bio import AlignIO

from primerforge.reference.reference import Reference
from primerforge.reference.coordinates import Coordinate
from primerforge.reference.exceptions import CoordinateSystemError


class CoordinateMapper:
    """
    Map coordinates between reference and alignment systems.
    """

    def __init__(
        self,
        reference: Reference,
    ):

        self.reference = reference

        if not reference.has_alignment:

            raise ValueError(
                "Reference has no alignment."
            )

        self.alignment = AlignIO.read(
            reference.alignment,
            "fasta",
        )

        self.reference_record = self.alignment[0]

        self.ref_to_aln = {}

        self.aln_to_ref = {}

        reference_position = 0

        for alignment_position, base in enumerate(
            str(self.reference_record.seq),
            start=1,
        ):

            if base != "-":

                reference_position += 1

                self.ref_to_aln[
                    reference_position
                ] = alignment_position

                self.aln_to_ref[
                    alignment_position
                ] = reference_position

    def reference_to_alignment(
        self,
        coordinate: Coordinate,
    ) -> Coordinate:

        if coordinate.system != "reference":

            raise CoordinateSystemError(
                self.reference.length,
                self.alignment.get_alignment_length(),
            )

        return Coordinate(

            start=self.ref_to_aln[
                coordinate.start
            ],

            end=self.ref_to_aln[
                coordinate.end
            ],

            system="alignment",

        )

    def alignment_to_reference(
        self,
        coordinate: Coordinate,
    ) -> Coordinate:

        if coordinate.system != "alignment":

            raise CoordinateSystemError(
                self.reference.length,
                self.alignment.get_alignment_length(),
            )

        return Coordinate(

            start=self.aln_to_ref[
                coordinate.start
            ],

            end=self.aln_to_ref[
                coordinate.end
            ],

            system="reference",

        )
