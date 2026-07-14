"""
Reference ↔ alignment coordinate conversion.

Primer locations are stored relative to the ungapped reference
sequence. This module converts those coordinates into alignment
columns so that primers can be evaluated correctly in gapped
multiple-sequence alignments.
"""

from pathlib import Path

from Bio import AlignIO


class CoordinateMapper:
    """
    Maps reference coordinates onto alignment coordinates.
    """

    def __init__(self, alignment: Path):

        self.alignment = AlignIO.read(
            alignment,
            "fasta",
        )

        self.reference = str(
            self.alignment[0].seq
        )

    def reference_to_alignment(
        self,
        position: int,
    ) -> int:
        """
        Convert a 0-based reference coordinate into an alignment
        column.
        """

        ungapped = 0

        for column, base in enumerate(self.reference):

            if base != "-":

                if ungapped == position:
                    return column

                ungapped += 1

        raise ValueError(
            f"Reference position {position} outside reference."
        )

    def alignment_to_reference(
        self,
        column: int,
    ) -> int:
        """
        Convert an alignment column back into a reference coordinate.
        """

        if column >= len(self.reference):
            raise ValueError("Column outside alignment.")

        ungapped = 0

        for i in range(column):

            if self.reference[i] != "-":
                ungapped += 1

        return ungapped
