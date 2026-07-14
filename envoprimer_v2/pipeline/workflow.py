"""
EnvoPrimer workflow.

High-level orchestration of the complete primer design pipeline.
"""

from __future__ import annotations

from pathlib import Path

from envoprimer_v2.pipeline.design import DesignPipeline
from envoprimer_v2.report.csv import CSVReport


class EnvoPrimerWorkflow:

    def __init__(self):

        self.pipeline = DesignPipeline()

        self.csv = CSVReport()

    def run(

        self,

        reference_sequence: str,

        target_alignment: Path,

        background_alignment: Path,

        background_sequences: dict,

        output_directory: Path,

    ):

        output_directory.mkdir(

            parents=True,

            exist_ok=True,

        )

        print()

        print("=" * 70)
        print("EnvoPrimer")
        print("=" * 70)

        pairs = self.pipeline.run(

            reference_sequence=reference_sequence,

            target_alignment_fasta=target_alignment,

            background_alignment_fasta=background_alignment,

            background_sequences=background_sequences,

        )

        if not pairs:

            print()

            print("No primer pairs passed all filters.")

            return []

        self.csv.write(

            pairs,

            output_directory / "primers.csv",

        )

        print()

        print("=" * 70)

        print("Finished")

        print("=" * 70)

        print()

        print(f"Accepted primer pairs : {len(pairs)}")

        print(f"Best score            : {pairs[0].final_score:.2f}")

        print()

        print(

            "CSV:",

            output_directory / "primers.csv",

        )

        return pairs
