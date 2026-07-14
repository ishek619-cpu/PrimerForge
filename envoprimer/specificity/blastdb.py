"""
BLAST database utilities.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


class BlastDatabase:
    """
    Create nucleotide BLAST databases.
    """

    def create(
        self,
        fasta: Path,
        output: Path,
    ) -> Path:

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        cmd = [

            "makeblastdb",

            "-in",
            str(fasta),

            "-dbtype",
            "nucl",

            "-out",
            str(output),

        ]

        subprocess.run(
            cmd,
            check=True,
        )

        return output

    def build(
        self,
        fasta: Path,
        output: str | Path,
    ) -> Path:
        """
        Backwards-compatible wrapper.
        """

        return self.create(
            fasta=fasta,
            output=Path(output),
        )

    def exists(
        self,
        database: Path,
    ) -> bool:

        extensions = (

            ".nhr",

            ".nin",

            ".nsq",

        )

        return all(

            (database.with_suffix(ext)).exists()

            for ext in extensions

        )
