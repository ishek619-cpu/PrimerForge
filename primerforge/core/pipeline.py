"""
Complete PrimerForge pipeline.
"""

from pathlib import Path

from primerforge.config.config import Config

from primerforge.io.sequence_manager import SequenceManager
from primerforge.io.merge import FASTAMerger

from primerforge.analysis.align import MAFFTAligner
from primerforge.analysis.snps import SNPFinder

from primerforge.primer.discovery import PrimerDiscovery

from primerforge.specificity.engine import SpecificityEngine

from primerforge.validation.validator import PrimerValidator

from primerforge.report.csv import CSVReport
from primerforge.report.json import JSONReport
from primerforge.report.html import HTMLReport
from primerforge.report.excel import ExcelReport


class Pipeline:
    """
    Complete PrimerForge workflow.
    """

    def __init__(self):

        self.sequence_manager = SequenceManager()

        self.merger = FASTAMerger()

        self.aligner = MAFFTAligner()

        self.snpfinder = SNPFinder()

        self.discovery = PrimerDiscovery()

        self.validator = PrimerValidator()

        self.csv = CSVReport()

        self.json = JSONReport()

        self.html = HTMLReport()

        self.excel = ExcelReport()

    def run(
        self,
        config_file: Path,
    ):

        config = Config(config_file)

        print("Downloading target sequences...")

        target = self.sequence_manager.get(
            config.organism,
            config.gene,
            config.download.get(
                "max_records",
                1000,
            ),
        )

        print("Downloading contrast taxa...")

        contrast = self.sequence_manager.get_contrast(
            config.contrast_taxa,
            config.gene,
            config.download.get(
                "max_records",
                1000,
            ),
        )

        results = Path("results")

        results.mkdir(
            exist_ok=True,
        )

        #
        # Merge FASTA files
        #

        target_merged = results / "target.fasta"

        self.merger.merge(
            [target],
            target_merged,
        )

        contrast_merged = results / "contrast.fasta"

        self.merger.merge(
            contrast,
            contrast_merged,
        )

        print()

        print(f"Target FASTA   : {target_merged}")

        print(f"Contrast FASTA : {contrast_merged}")

        alignment = results / "alignment.fasta"

        print()

        print("Running MAFFT...")

        self.aligner.align(
            target_merged,
            alignment,
        )

        print("Finding SNPs...")

        snps = self.snpfinder.find(
            alignment,
        )

        #
        # Build BLAST database
        #

        print("Building BLAST database...")

        blast_db = "results/blast/primerforge"

        self.discovery.blastdb.build(
            alignment,
            blast_db,
        )

        #
        # Create specificity engine
        #

        print("Creating specificity engine...")

        engine = SpecificityEngine(
            database=blast_db,
            target_species=config.organism,
        )

        #
        # Discover primers
        #

        print("Designing primers...")

        pairs = self.discovery.discover(
            target_merged,
            specificity_engine=engine,
        )

        validated = []

        print("Validating primers...")

        for pair in pairs:

            self.validator.validate(
                pair,
                alignment,
                snps,
            )

            validated.append(
                pair,
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

        print(f"Species      : {config.organism}")

        print(f"Marker       : {config.gene}")

        print(f"Target FASTA : {target_merged}")

        print(f"Contrast     : {len(contrast)} datasets")

        print(f"Alignment    : {alignment}")

        print(f"Primer pairs : {len(validated)}")

        if validated:

            print(
                f"Best score   : {validated[0].score:.2f}"
            )

        print()

        print("Reports")

        print("-------")

        print("CSV   :", results / "primers.csv")

        print("JSON  :", results / "primers.json")

        print("HTML  :", results / "index.html")

        print("Excel :", results / "primers.xlsx")

        print()
