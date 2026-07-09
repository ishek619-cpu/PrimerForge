"""
PrimerForge complete pipeline.
"""

from pathlib import Path

from primerforge.config.config import Config
from primerforge.primer.discovery import PrimerDiscovery
from primerforge.analysis.snps import SNPFinder
from primerforge.validation.validator import PrimerValidator
from primerforge.report.csv import CSVReport


class Pipeline:

    def __init__(
        self,
        config_file,
    ):

        self.config = Config(
            Path(config_file),
        )

        self.discovery = PrimerDiscovery()

        self.validator = PrimerValidator()

        self.snpfinder = SNPFinder()

        self.report = CSVReport()

    def run(self):

        reference = Path(
            "data/genes/NC_013663_CYTB.fasta"
        )

        alignment = Path(
            "data/alignments/alignment.fasta"
        )

        pairs = self.discovery.discover(
            reference,
        )

        snps = self.snpfinder.find(
            alignment,
        )

        validated = []

        for pair in pairs:

            result = self.validator.validate(
                pair,
                alignment,
                snps,
            )

            pair.score = result["final_score"]

            validated.append(pair)

        validated.sort(
            key=lambda x: x.score,
            reverse=True,
        )

        self.report.write(
            validated,
            Path(
                "results/primers.csv"
            ),
        )

        print()

        print("=" * 60)
        print("PrimerForge completed successfully")
        print("=" * 60)
        print()

        print(f"Primer pairs : {len(validated)}")
        print(f"Best score   : {validated[0].score:.2f}")
        print("CSV report   : results/primers.csv")
        print()
