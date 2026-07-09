"""
Complete PrimerForge pipeline.
"""

from pathlib import Path

from primerforge.config.config import Config
from primerforge.primer.discovery import PrimerDiscovery
from primerforge.analysis.snps import SNPFinder
from primerforge.validation.validator import PrimerValidator

from primerforge.report.csv import CSVReport
from primerforge.report.json import JSONReport
from primerforge.report.html import HTMLReport
from primerforge.report.excel import ExcelReport


class Pipeline:
    """
    Main PrimerForge pipeline.
    """

    def __init__(self):

        self.discovery = PrimerDiscovery()

        self.snpfinder = SNPFinder()

        self.validator = PrimerValidator()

        self.csv = CSVReport()

        self.json = JSONReport()

        self.html = HTMLReport()

        self.excel = ExcelReport()

    def run(
        self,
        config_file: Path,
    ):

        config = Config(
            config_file,
        )

        gene = Path(
            "data/genes/NC_013663_CYTB.fasta"
        )

        alignment = Path(
            "data/alignments/alignment.fasta"
        )

        pairs = self.discovery.discover(
            gene,
        )

        snps = self.snpfinder.find(
            alignment,
        )

        validated = []

        for pair in pairs:

            self.validator.validate(
                pair,
                alignment,
                snps,
            )

            validated.append(
                pair,
            )

        results = Path(
            "results"
        )

        self.csv.write(
            validated,
            results / "primers.csv",
        )

        self.json.write(
            validated,
            results / "primers.json",
        )

        self.html.write(
            validated,
            results / "index.html",
        )

        self.excel.write(
            validated,
            results / "primers.xlsx",
        )

        print()

        print("=" * 60)

        print("PrimerForge completed successfully")

        print("=" * 60)

        print()

        print(
            f"Primer pairs : {len(validated)}"
        )

        if validated:

            print(
                f"Best score   : {validated[0].score:.2f}"
            )

        print()

        print("Reports generated")

        print("-----------------")

        print("CSV   : results/primers.csv")

        print("JSON  : results/primers.json")

        print("HTML  : results/index.html")

        print("Excel : results/primers.xlsx")

        print()
