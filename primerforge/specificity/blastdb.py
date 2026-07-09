"""
Create local BLAST databases.
"""

import subprocess
from pathlib import Path


class BlastDatabase:
    """
    Build a nucleotide BLAST database.
    """

    def build(
        self,
        fasta: Path,
        database: str,
    ):

        cmd = [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "nucl",
            "-out",
            database,
        ]

        subprocess.run(
            cmd,
            check=True,
        )

        return database
