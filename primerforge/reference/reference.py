"""
Reference sequence model.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass(slots=True)
class Reference:

    fasta: Path

    alignment: Path | None = None

    annotation: Path | None = None

    accession: str = ""

    description: str = ""

    sequence: str = ""

    length: int = 0

    def __post_init__(self):

        record = next(
            SeqIO.parse(
                self.fasta,
                "fasta",
            )
        )

        self.accession = record.id

        self.description = record.description

        self.sequence = str(
            record.seq
        )

        self.length = len(
            self.sequence
        )

    @property
    def has_alignment(self):

        return self.alignment is not None

    @property
    def has_annotation(self):

        return self.annotation is not None
