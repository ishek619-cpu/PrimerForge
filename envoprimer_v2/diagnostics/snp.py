"""
Diagnostic SNP detection.

A diagnostic SNP is a position where all target sequences share one
allele and every background sequence possesses a different allele.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer_v2.models.alignment import Alignment


@dataclass(slots=True)
class DiagnosticSNP:

    alignment_position: int

    target_base: str

    background_bases: set[str]

    target_support: float

    background_support: float

    score: float


class DiagnosticSNPFinder:

    def __init__(

        self,

        minimum_target_support: float = 1.0,

        maximum_background_support: float = 0.0,

    ):

        self.minimum_target_support = minimum_target_support

        self.maximum_background_support = maximum_background_support

    def find(

        self,

        target: Alignment,

        background: Alignment,

    ) -> list[DiagnosticSNP]:

        if target.length != background.length:

            raise ValueError(

                "Target and background alignments have different lengths."

            )

        diagnostics = []

        for column in range(target.length):

            target_column = [

                b.upper()

                for b in target.column(column)

                if b != "-"

            ]

            background_column = [

                b.upper()

                for b in background.column(column)

                if b != "-"

            ]

            if not target_column:

                continue

            if not background_column:

                continue

            target_base = max(

                set(target_column),

                key=target_column.count,

            )

            target_support = (

                target_column.count(target_base)

                / len(target_column)

            )

            if target_support < self.minimum_target_support:

                continue

            background_bases = set(

                background_column

            )

            if target_base in background_bases:

                continue

            background_support = 0.0

            score = (

                target_support

                * 100.0

            )

            diagnostics.append(

                DiagnosticSNP(

                    alignment_position=column,

                    target_base=target_base,

                    background_bases=background_bases,

                    target_support=round(

                        target_support,

                        4,

                    ),

                    background_support=background_support,

                    score=round(

                        score,

                        2,

                    ),

                )

            )

        diagnostics.sort(

            key=lambda x: (

                -x.score,

                x.alignment_position,

            )

        )

        return diagnostics
