"""
Local BLAST wrapper.
"""

import subprocess
from pathlib import Path


class BlastRunner:
    """
    Execute blastn-short searches.
    """

    def __init__(self, threads: int = 4):

        self.threads = threads

    def search(
        self,
        query: Path,
        database: str,
        output: Path,
    ) -> Path:

        cmd = [

            "blastn",

            "-task",
            "blastn-short",

            "-query",
            str(query),

            "-db",
            database,

            "-outfmt",
            "6",

            "-num_threads",
            str(self.threads),

            "-out",
            str(output),
        ]

        subprocess.run(
            cmd,
            check=True,
        )

        return output
