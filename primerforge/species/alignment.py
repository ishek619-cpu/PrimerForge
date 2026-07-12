"""
Species alignment workflow.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from primerforge.species.aligner import MAFFTAligner


@dataclass(slots=True)
class AlignmentResult:
    """
    Result of species alignments.
    """

    target_alignment: Path

    background_alignment: Path


class AlignmentBuilder:
    """
    Build MAFFT alignments for target and background datasets.
    """

    def __init__(
        self,
        executable: str = "mafft",
    ):

        self.aligner = MAFFTAligner(
            executable=executable,
        )

    def build(
        self,
        target_fasta: Path,
        background_fasta: Path,
        output_directory: Path,
    ) -> AlignmentResult:

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        target_alignment = (
            output_directory
            / "target_alignment.fasta"
        )

        background_alignment = (
            output_directory
            / "background_alignment.fasta"
        )

        print("Aligning target sequences...")

        self.aligner.align(
            target_fasta,
            target_alignment,
        )

        print("Aligning background sequences...")

        self.aligner.align(
            background_fasta,
            background_alignment,
        )

        return AlignmentResult(

            target_alignment=target_alignment,

            background_alignment=background_alignment,

        )
