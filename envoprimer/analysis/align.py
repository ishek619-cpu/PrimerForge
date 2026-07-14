"""
Multiple sequence alignment using MAFFT.
"""

import subprocess
from pathlib import Path


class MAFFTAligner:
    """
    Run MAFFT and produce an aligned FASTA.
    """

    def align(
        self,
        sequences: Path,
        output: Path,
        threads: int = -1,
    ) -> Path:

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with open(
            output,
            "w",
            encoding="utf-8",
        ) as handle:

            subprocess.run(
                [
                    "mafft",
                    "--thread",
                    str(threads),
                    "--auto",
                    str(sequences),
                ],
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

        return output
