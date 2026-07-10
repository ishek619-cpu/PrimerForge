"""
Complete PrimerForge pipeline.
"""

from pathlib import Path

from primerforge.config.config import Config
from primerforge.io.sequence_manager import SequenceManager
from primerforge.analysis.align import MAFFTAligner
from primerforge.analysis.snps import SNPFinder
from primerforge.primer.discovery import PrimerDiscovery
from primerforge.validation.validator import PrimerValidator
from primerforge.report.csv import CSVReport
from primerforge.report.json import JSONReport
from primerforge.report.html import HTMLReport
from primerforge.report.excel import ExcelReport


class Pipeline:

    def __init__(self):

        self.sequence_manager = SequenceManager()
        self.aligner = MAFFTAligner()
        self.snpfinder = SNPFinder()
        self.discovery = PrimerDiscovery()
        self.validator = PrimerValidator()

        self.csv = CSVReport()
        self.json = JSONReport()
        self.html = HTMLReport()
        self.excel = ExcelReport()

    def run(self, config_file: Path):

        print("STEP 1 - Loading config")

        config = Config(config_file)

        print("STEP 2 - Download/cache sequences")

        gene = self.sequence_manager.get(
            config.organism,
            config.gene,
        )

        print("STEP 3 - Alignment")

        alignment = Path("results/alignment.fasta")

        self.aligner.align(
            gene,
            alignment,
        )

        print("STEP 4 - SNP discovery")

        snps = self.snpfinder.find(
            alignment,
        )

        print("STEP 5 - Primer discovery")

        pairs = self.discovery.discover(
            gene,
        )

        print("STEP 6 - Validation")

        validated = []

        for pair in pairs:

            self.validator.validate(
                pair,
                alignment,
                snps,
            )

            validated.append(pair)

        print("STEP 7 - Reports")

        results = Path("results")
        results.mkdir(exist_ok=True)

        self.csv.write(validated, results / "primers.csv")
        self.json.write(validated, results / "primers.json")
        self.html.write(validated, results / "index.html")
        self.excel.write(validated, results / "primers.xlsx")

        print("DONE")
