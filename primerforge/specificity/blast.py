"""
Local BLAST wrapper.
"""

import subprocess
from pathlib import Path


class BlastRunner:
    """
    Execute blastn-short searches.
    """

    def __init__(
        self,
        threads: int = 4,
    ):

        self.threads = threads

    def search(
        self,
        query: Path,
        database: str,
        output: Path,
    ) -> Path:

        outfmt = (
            "6 "
            "qseqid "
            "sseqid "
            "pident "
            "length "
            "mismatch "
            "gapopen "
            "qstart "
            "qend "
            "sstart "
            "send "
            "qseq "
            "sseq "
            "evalue "
            "bitscore "
            "sstrand"
        )

        cmd = [

            "blastn",

            "-task",
            "blastn-short",

            "-query",
            str(query),

            "-db",
            database,

            "-outfmt",
            outfmt,

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
