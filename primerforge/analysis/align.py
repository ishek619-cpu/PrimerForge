"""
PrimerForge Alignment Engine

Automatically combines FASTA files and aligns them with MUSCLE.
"""

from pathlib import Path
import subprocess

from primerforge.core.logger import get_logger


class AlignmentEngine:

    def __init__(self):

        self.logger = get_logger(__name__)

    def combine_fastas(
        self,
        fasta_dir: Path,
        output_fasta: Path,
    ) -> Path:

        fasta_files = sorted(
            fasta_dir.glob("*.fasta")
        )

        if len(fasta_files) == 0:
            raise FileNotFoundError(
                "No FASTA files found."
            )

        output_fasta.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with open(output_fasta, "w") as outfile:

            for fasta in fasta_files:

                self.logger.info(
                    f"Adding {fasta.name}"
                )

                with open(fasta) as infile:

                    outfile.write(infile.read())

                    if not infile.read().endswith("\n"):
                        outfile.write("\n")

        self.logger.info(
            f"Combined {len(fasta_files)} FASTA files"
        )

        return output_fasta

    def run_muscle(
        self,
        input_fasta: Path,
        output_fasta: Path,
    ) -> Path:

        self.logger.info(
            f"Running MUSCLE"
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
