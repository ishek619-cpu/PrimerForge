"""
Species-specific PrimerForge pipeline.
"""

from __future__ import annotations

from pathlib import Path

from primerforge.config.config import Config
from primerforge.report.csv import CSVReport
from primerforge.report.excel import ExcelReport
from primerforge.report.html import HTMLReport
from primerforge.report.json import JSONReport
from primerforge.species.designer import SpeciesPrimerDesigner
from primerforge.species.workflow import SpeciesWorkflow


class SpeciesPipeline:
    """
    Complete species-specific primer design pipeline.
    """

    def __init__(self):

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

        results = Path("results")

        results.mkdir(
            exist_ok=True,
        )

        print()

        print("=" * 60)
        print("PrimerForge Species Pipeline")
        print("=" * 60)

        workflow = SpeciesWorkflow(

            email=config.ncbi_email,

            api_key=config.ncbi_api_key,

        )

        workflow_result = workflow.run(

            species=config.organism,

            marker=config.gene,

            output_directory=results,

        )

        print()

        print(
            f"Diagnostic SNPs: {len(workflow_result.diagnostic_sites)}"
        )

        print(
            f"Diagnostic windows: {len(workflow_result.diagnostic_windows)}"
        )

        designer = SpeciesPrimerDesigner()

        print()

        print("Designing primers...")

        pairs = designer.design(
            workflow_result,
        )

        self.csv.write(
            pairs,
            results / "primers.csv",
        )

        self.json.write(
            pairs,
            results / "primers.json",
        )

        self.html.write(
            pairs,
            results / "index.html",
        )

        self.excel.write(
            pairs,
            results / "primers.xlsx",
        )

        print()

        print("=" * 60)
        print("Species pipeline completed")
        print("=" * 60)

        print()

        print(
            f"Primer pairs: {len(pairs)}"
        )

        if pairs:

            print(
                f"Best score: {pairs[0].score:.2f}"
            )
