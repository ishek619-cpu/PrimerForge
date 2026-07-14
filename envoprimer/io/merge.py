"""
Merge FASTA files.
"""

from pathlib import Path
from Bio import SeqIO


class FASTAMerger:
    """
    Merge multiple FASTA files into one.
    """

    def merge(
        self,
        files: list[Path],
        output: Path,
    ) -> Path:

        records = []

        seen = set()

        for fasta in files:

            for record in SeqIO.parse(
                fasta,
                "fasta",
            ):

                sequence = str(record.seq)

                if sequence in seen:
                    continue

                seen.add(sequence)

                records.append(record)

        SeqIO.write(
            records,
            output,
            "fasta",
        )

        return output
