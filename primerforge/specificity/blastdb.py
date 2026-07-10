"""
BLAST database utilities.
"""

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
