"""
Population conservation analysis.

Calculates conservation at every alignment column and
provides window-level conservation scores.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer_v2.models.alignment import Alignment


@dataclass(slots=True)
class ConservedSite:

    alignment_position: int

    consensus: str

    conservation: float

    depth: int


class ConservationAnalyzer:

    def analyse(

        self,

        alignment: Alignment,

    ) -> list[ConservedSite]:

        results = []

        for column in range(alignment.length):

            bases = [

                b.upper()

                for b in alignment.column(column)

                if b != "-"

            ]

            if not bases:

                continue

            consensus = max(

                set(bases),

                key=bases.count,

            )

            conservation = (

                bases.count(consensus)

                / len(bases)

            ) * 100.0

            results.append(

                ConservedSite(

                    alignment_position=column,

                    consensus=consensus,

                    conservation=round(

                        conservation,

                        2,

                    ),

                    depth=len(bases),

                )

            )

        return results

    def window_score(

        self,

        sites: list[ConservedSite],

        start: int,

        end: int,

    ) -> float:

        window = [

            site.conservation

            for site in sites

            if start <= site.alignment_position <= end

        ]

        if not window:

            return 0.0

        return round(

            sum(window)

            / len(window),

            2,

        )
