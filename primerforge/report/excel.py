"""
Excel report generator.
"""

from pathlib import Path

from openpyxl import Workbook

from primerforge.models.pair import PrimerPair


class ExcelReport:
    """
    Generate Excel report.
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

        workbook = Workbook()

        sheet = workbook.active

        sheet.title = "Primers"

        sheet.append([
            "Forward",
            "Reverse",
            "Forward Tm",
            "Reverse Tm",
            "Forward GC",
            "Reverse GC",
            "Product Size",
            "Score",
        ])

        for pair in pairs:

            sheet.append([

                pair.forward.sequence,

                pair.reverse.sequence,

                round(pair.forward.tm,2),

                round(pair.reverse.tm,2),

                pair.forward.gc,

                pair.reverse.gc,

                pair.product_size,

                pair.score,

            ])

        workbook.save(output)
