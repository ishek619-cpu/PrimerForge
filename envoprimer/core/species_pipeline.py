"""
Species-specific EnvoPrimer pipeline.
"""

from __future__ import annotations

import json

from pathlib import Path

from envoprimer.config.config import Config

from envoprimer.report.csv import CSVReport
from envoprimer.report.dataset_summary import DatasetSummaryReport
from envoprimer.report.excel import ExcelReport
from envoprimer.report.html import HTMLReport
from envoprimer.report.json import JSONReport

from envoprimer.species.designer import SpeciesPrimerDesigner
from envoprimer.species.workflow import SpeciesWorkflow


class SpeciesPipeline:
    """
    Complete species-specific primer design pipeline.
    """

    def __init__(self):

        self.csv = CSVReport()

        self.json = JSONReport()

        self.html = HTMLReport()

        self.excel = ExcelReport()

        self.dataset = DatasetSummaryReport()

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
        print("EnvoPrimer Species Pipeline")
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

        #
        # EnvoPrimer Engine
        #

        designer = SpeciesPrimerDesigner()

        print()

        print("=" * 60)
        print("Running EnvoPrimer Engine")
        print("=" * 60)

        print()

        pairs = designer.design(
            workflow_result,
        )

        print()

        print(
            f"Accepted primer pairs : {len(pairs)}"
        )

        #
        # Primer reports
        #

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

        #
        # Dataset reports
        #

        metadata = (
            workflow_result.dataset.target_metadata
            + workflow_result.dataset.background_metadata
        )

        self.dataset.write_csv(
            metadata,
            results / "dataset_summary.csv",
        )

        self.dataset.write_excel(
            metadata,
            results / "dataset_summary.xlsx",
        )

        #
        # Accession list
        #

        accession_file = (
            results / "accessions.txt"
        )

        with accession_file.open(
            "w",
            encoding="utf-8",
        ) as handle:

            handle.write(
                "TARGET SEQUENCES\n"
            )

            handle.write(
                "================\n\n"
            )

            for record in workflow_result.dataset.target_metadata:

                handle.write(
                    f"{record.accession}\t{record.species}\n"
                )

            handle.write(
                "\nBACKGROUND SEQUENCES\n"
            )

            handle.write(
                "====================\n\n"
            )

            for record in workflow_result.dataset.background_metadata:

                handle.write(
                    f"{record.accession}\t{record.species}\n"
                )

        #
        # Workflow metadata
        #

        workflow_metadata = {

            "species": workflow_result.dataset.species,

            "marker": workflow_result.dataset.marker,

            "target_sequences": len(
                workflow_result.dataset.target_metadata
            ),

            "background_sequences": len(
                workflow_result.dataset.background_metadata
            ),

            "diagnostic_snps": len(
                workflow_result.diagnostic_sites
            ),

            "diagnostic_windows": len(
                workflow_result.diagnostic_windows
            ),

            "primer_pairs": len(
                pairs
            ),

        }

        with (
            results
            / "workflow_metadata.json"
        ).open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                workflow_metadata,
                handle,
                indent=4,
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
                f"Best score: {pairs[0].final_score:.2f}"
            )

        print()

        print("Reports")
        print("-------")

        print(
            "CSV:",
            results / "primers.csv",
        )

        print(
            "JSON:",
            results / "primers.json",
        )

        print(
            "HTML:",
            results / "index.html",
        )

        print(
            "Excel:",
            results / "primers.xlsx",
        )

        print(
            "Dataset CSV:",
            results / "dataset_summary.csv",
        )

        print(
            "Dataset Excel:",
            results / "dataset_summary.xlsx",
        )

        print(
            "Accessions:",
            results / "accessions.txt",
        )

        print(
            "Workflow:",
            results / "workflow_metadata.json",
        )

        print()
