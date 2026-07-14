"""
Sequence filtering utilities.
"""

from pathlib import Path

from Bio import SeqIO


class SequenceFilter:
    """
    Filter FASTA sequences by minimum length.
    """

    def filter_by_length(
        self,
        input_fasta: Path,
        output_fasta: Path,
        minimum_length: int = 1000,
    ) -> Path:

        records = [
            record
            for record in SeqIO.parse(input_fasta, "fasta")
            if len(record.seq) >= minimum_length
        ]

        output_fasta.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        SeqIO.write(
            records,
            output_fasta,
            "fasta",
        )

        print(
            f"Retained {len(records)} sequences >= {minimum_length} bp"
        )

        return output_fasta
