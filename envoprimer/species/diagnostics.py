"""
Diagnostic SNP discovery.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio import AlignIO


@dataclass(slots=True)
class DiagnosticSite:

    position: int

    target_base: str

    background_bases: set[str]


class DiagnosticFinder:

    def find(
        self,
        target_alignment,
        background_alignment,
    ) -> list[DiagnosticSite]:

        target = AlignIO.read(
            target_alignment,
            "fasta",
        )

        background = AlignIO.read(
            background_alignment,
            "fasta",
        )

        reference = str(
            target[0].seq,
        )

        #
        # Alignment column -> reference coordinate
        #
        column_to_reference = {}

        reference_position = 0

        for column, base in enumerate(reference):

            if base != "-":

                reference_position += 1

                column_to_reference[column] = reference_position - 1

        length = min(

            target.get_alignment_length(),

            background.get_alignment_length(),

        )

        sites = []

        for column in range(length):

            #
            # Ignore columns that are gaps in the reference
            #
            if column not in column_to_reference:

                continue

            target_bases = {

                b

                for b in target[:, column]

                if b != "-"

            }

            background_bases = {

                b

                for b in background[:, column]

                if b != "-"

            }

            if len(target_bases) != 1:

                continue

            base = next(
                iter(target_bases),
            )

            if base in background_bases:

                continue

            sites.append(

                DiagnosticSite(

                    position=column_to_reference[column],

                    target_base=base,

                    background_bases=background_bases,

                )

            )

        return sites
