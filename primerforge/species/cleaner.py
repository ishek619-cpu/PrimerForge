"""
Sequence cleaning utilities.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO


class SequenceCleaner:

    def clean(
        self,
        input_fasta: Path,
        output_fasta: Path,
        min_length: int = 100,
        max_n: float = 0.05,
    ) -> int:

        output_fasta.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        seen = set()

        kept = []

        for record in SeqIO.parse(
            input_fasta,
            "fasta",
        ):

            sequence = str(
                record.seq,
            ).upper()

            if len(sequence) < min_length:

                continue

            if (
                sequence.count("N")
                / len(sequence)
            ) > max_n:

                continue

            if sequence in seen:

                continue

            seen.add(sequence)

            kept.append(record)

        SeqIO.write(
            kept,
            output_fasta,
            "fasta",
        )

        return len(
            kept,
        )
