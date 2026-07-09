"""
CSV report generator.
"""

import csv
from pathlib import Path

from primerforge.models.pair import PrimerPair


class CSVReport:

    """
    Export primer pairs to CSV.
    """

    def write(
        self,
        pairs: list[PrimerPair],
        output: Path,
    ):

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with open(
            output,
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:

            writer = csv.writer(handle)

            writer.writerow(
                [
                    "Rank",
                    "Score",
                    "Forward",
                    "Reverse",
                    "Tm(F)",
                    "Tm(R)",
                    "GC(F)",
                    "GC(R)",
                    "Product",
                ]
            )

            for rank, pair in enumerate(
                pairs,
                start=1,
            ):

                writer.writerow(
                    [
                        rank,
                        pair.score,
                        pair.forward.sequence,
                        pair.reverse.sequence,
                        round(pair.forward.tm, 2),
                        round(pair.reverse.tm, 2),
                        round(pair.forward.gc, 2),
                        round(pair.reverse.gc, 2),
                        pair.product_size,
                    ]
                )

        return output
