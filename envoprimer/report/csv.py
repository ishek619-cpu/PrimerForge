"""
CSV report generator.
"""

import csv
from pathlib import Path


class CSVReport:
    """
    Export primer pairs to CSV.
    """

    def write(
        self,
        pairs,
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
                    "Final Score",
                    "Forward Primer",
                    "Reverse Primer",
                    "Tm(F)",
                    "Tm(R)",
                    "GC(F)",
                    "GC(R)",
                    "Product Size",
                    "Population Score",
                    "Specificity Score",
                    "Thermo Score",
                    "3Prime Score",
                    "BLAST Score",
                    "PCR Score",
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

                        round(pair.forward.tm, 2),
                        round(pair.reverse.tm, 2),

                        round(pair.forward.gc, 2),
                        round(pair.reverse.gc, 2),

                        pair.product_size,

                        round(pair.population_score, 2),
                        round(pair.specificity_score, 2),
                        round(pair.thermo_score, 2),
                        round(pair.three_prime_score, 2),
                        round(pair.blast_score, 2),
                        round(pair.pcr_score, 2),
                    ]
                )

        return output
