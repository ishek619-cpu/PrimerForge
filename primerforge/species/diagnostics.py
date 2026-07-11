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

        length = min(
            target.get_alignment_length(),
            background.get_alignment_length(),
        )

        sites = []

        for column in range(length):

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

            if base not in background_bases:

                sites.append(

                    DiagnosticSite(

                        position=column,

                        target_base=base,

                        background_bases=background_bases,

                    )

                )

        return sites
