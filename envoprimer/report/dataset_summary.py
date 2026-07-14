"""
Dataset summary report.
"""

from __future__ import annotations

import csv
from pathlib import Path

from openpyxl import Workbook

from envoprimer.species.metadata import SequenceMetadata


class DatasetSummaryReport:

    def write_csv(
        self,
        metadata: list[SequenceMetadata],
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
                    "Role",
                    "Species",
                    "Accession",
                    "Marker",
                    "Length",
                    "Definition",
                    "Source",
                ]
            )

            for record in metadata:

                writer.writerow(
                    [
                        record.role,
                        record.species,
                        record.accession,
                        record.marker,
                        record.length,
                        record.definition,
                        record.source,
                    ]
                )

    def write_excel(
        self,
        metadata: list[SequenceMetadata],
        output: Path,
    ):

        wb = Workbook()

        ws = wb.active

        ws.title = "Dataset"

        ws.append(
            [
                "Role",
                "Species",
                "Accession",
                "Marker",
                "Length",
                "Definition",
                "Source",
            ]
        )

        for record in metadata:

            ws.append(
                [
                    record.role,
                    record.species,
                    record.accession,
                    record.marker,
                    record.length,
                    record.definition,
                    record.source,
                ]
            )

        wb.save(output)
