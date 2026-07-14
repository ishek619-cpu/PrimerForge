"""
Species alignment workflow with persistent cache.
"""

from __future__ import annotations

import shutil

from dataclasses import dataclass
from pathlib import Path

from envoprimer.cache.engine import CacheEngine
from envoprimer.species.aligner import MAFFTAligner


@dataclass(slots=True)
class AlignmentResult:
    """
    Result of species alignments.
    """

    target_alignment: Path
    background_alignment: Path


class AlignmentBuilder:
    """
    Build MAFFT alignments with persistent caching.
    """

    def __init__(
        self,
        executable: str = "mafft",
    ):

        self.aligner = MAFFTAligner(
            executable=executable,
        )

        self.cache = CacheEngine()

    def build(
        self,
        target_fasta: Path,
        background_fasta: Path,
        output_directory: Path,
        marker: str,
    ) -> AlignmentResult:

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        species = target_fasta.parent.name.replace(
            "_",
            " ",
        )

        marker = marker.upper()

        cache_dir = self.cache.alignments(
            species,
            marker,
        )

        cached_target = (
            cache_dir
            / "target_alignment.fasta"
        )

        cached_background = (
            cache_dir
            / "background_alignment.fasta"
        )

        target_alignment = (
            output_directory
            / "target_alignment.fasta"
        )

        background_alignment = (
            output_directory
            / "background_alignment.fasta"
        )

        ####################################################
        # TARGET
        ####################################################

        if cached_target.exists():

            print()
            print("=" * 60)
            print("USING CACHED TARGET ALIGNMENT")
            print("=" * 60)

            shutil.copy2(
                cached_target,
                target_alignment,
            )

        else:

            print("Aligning target sequences...")

            self.aligner.align(
                target_fasta,
                target_alignment,
            )

            shutil.copy2(
                target_alignment,
                cached_target,
            )

        ####################################################
        # BACKGROUND
        ####################################################

        if cached_background.exists():

            print()
            print("=" * 60)
            print("USING CACHED BACKGROUND ALIGNMENT")
            print("=" * 60)

            shutil.copy2(
                cached_background,
                background_alignment,
            )

        else:

            print("Aligning background sequences...")

            self.aligner.align(
                background_fasta,
                background_alignment,
            )

            shutil.copy2(
                background_alignment,
                cached_background,
            )

        return AlignmentResult(

            target_alignment=target_alignment,

            background_alignment=background_alignment,

        )
