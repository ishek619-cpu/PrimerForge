"""
PrimerForge Alignment Engine

Runs MUSCLE on FASTA files.
"""

from pathlib import Path
import subprocess

from primerforge.core.logger import get_logger


class AlignmentEngine:
    """
    Perform multiple sequence alignments using MUSCLE.
    """

    def __init__(self):

        self.logger = get_logger(__name__)

    def run_muscle(
        self,
        input_fasta: Path,
        output_fasta: Path,
    ) -> Path:
        """
        Run MUSCLE alignment.
        """

        self.logger.info(
            f"Running MUSCLE on {input_fasta}"
        )

        output_fasta.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        command = [
            "muscle",
            "-align",
            str(input_fasta),
            "-output",
            str(output_fasta),
        ]

        subprocess.run(
            command,
            check=True,
        )

        self.logger.info(
            f"Alignment saved to {output_fasta}"
        )

        return output_fasta

