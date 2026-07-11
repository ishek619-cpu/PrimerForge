"""
MAFFT alignment wrapper.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


class MAFFTAligner:

    def __init__(
        self,
        executable: str = "mafft",
    ):

        self.executable = executable

    def align(
        self,
        input_fasta: Path,
        output_fasta: Path,
        threads: int = -1,
    ) -> Path:

        output_fasta.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with output_fasta.open(
            "w",
            encoding="utf-8",
        ) as out:

            subprocess.run(

                [
                    self.executable,
                    "--thread",
                    str(threads),
                    "--auto",
                    str(input_fasta),
                ],

                stdout=out,

                check=True,

            )

        return output_fasta
