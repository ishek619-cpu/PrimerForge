"""
CSV export for EnvoPrimer.
"""

from __future__ import annotations

import csv
from pathlib import Path


class CSVReport:

    def write(

        self,

        pairs,

        output: Path,

    ):

        output.parent.mkdir(

            parents=True,

            exist_ok=True,

        )

        with output.open(

            "w",

            newline="",

            encoding="utf-8",

        ) as handle:

            writer = csv.writer(handle)

            writer.writerow(

                [

                    "Rank",

                    "Final Score",

                    "Forward Primer",

                    "Reverse Primer",

                    "Product Size",

                    "Forward Tm",

                    "Reverse Tm",

                    "Forward GC",

                    "Reverse GC",

                    "Primer3 Penalty",

                ]

            )

            for rank, pair in enumerate(

                pairs,

                start=1,

            ):

                writer.writerow(

                    [

                        rank,

                        round(pair.final_score, 2),

                        pair.forward.sequence,

                        pair.reverse.sequence,

                        pair.product_size,

                        round(pair.forward.tm, 2),

                        round(pair.reverse.tm, 2),

                        round(pair.forward.gc, 2),

                        round(pair.reverse.gc, 2),

                        round(pair.primer3_penalty, 3),

                    ]

                )

        return output
