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

                    "Product(bp)",

                    "Population Conservation",
                    "Population Coverage",

                    "Specificity Score",
                    "Specificity Passed",

                    "Thermo",
                    "Conservation",
                    "SNP",
                    "Product",
                    "Tm Balance",
                    "GC Balance",
                    "Multiplex",
                ]
            )

            for rank, pair in enumerate(
                pairs,
                start=1,
            ):

                breakdown = pair.breakdown

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

                        round(
                            pair.population_conservation,
                            2,
                        ),

                        round(
                            pair.population_coverage,
                            2,
                        ),

                        round(
                            pair.specificity_score,
                            2,
                        ),

                        pair.passed_specificity,

                        breakdown.get(
                            "thermo",
                            "",
                        ),

                        breakdown.get(
                            "conservation",
                            "",
                        ),

                        breakdown.get(
                            "snp",
                            "",
                        ),

                        breakdown.get(
                            "product",
                            "",
                        ),

                        breakdown.get(
                            "tm_balance",
                            "",
                        ),

                        breakdown.get(
                            "gc_balance",
                            "",
                        ),

                        breakdown.get(
                            "multiplex",
                            "",
                        ),
                    ]
                )

        return output
