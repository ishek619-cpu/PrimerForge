
"""
Sequence cleaning utilities.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO


class SequenceCleaner:
    """
    Clean downloaded marker datasets before alignment.
    """

    BAD_KEYWORDS = (
        "whole genome",
        "chromosome",
        "scaffold",
        "contig",
        "wgs",
        "genome assembly",
        "metagenome",
        "environmental",
        "uncultured",
        "synthetic",
        "vector",
        "plasmid",
        "clone",
    )

    def clean(
        self,
        input_fasta: Path,
        output_fasta: Path,
        min_length: int = 100,
        max_length: int = 3000,
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

            description = record.description.lower()

            #
            # Remove unwanted records
            #
            if any(
                keyword in description
                for keyword in self.BAD_KEYWORDS
            ):
                continue

            sequence = str(
                record.seq,
            ).upper()

            #
            # Length filter
            #
            if len(sequence) < min_length:
                continue

            if len(sequence) > max_length:
                continue

            #
            # Ambiguous bases
            #
            if (
                sequence.count("N")
                / len(sequence)
            ) > max_n:
                continue

            #
            # Remove duplicate sequences
            #
            if sequence in seen:
                continue

            seen.add(sequence)

            kept.append(record)

        SeqIO.write(
            kept,
            output_fasta,
            "fasta",
        )

        return len(kept)
