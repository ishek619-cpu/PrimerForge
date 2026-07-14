"""
Sequence model for EnvoPrimer.

Each sequence stores both aligned and ungapped forms together
with fast coordinate mappings.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class Sequence:

    id: str

    species: str

    accession: str

    aligned: str

    ungapped: str = ""

    alignment_to_sequence: list[int] = field(default_factory=list)

    sequence_to_alignment: list[int] = field(default_factory=list)

    def __post_init__(self):

        self.ungapped = self.aligned.replace("-", "")

        self.alignment_to_sequence = []

        self.sequence_to_alignment = []

        sequence_position = 0

        for alignment_position, base in enumerate(self.aligned):

            if base == "-":

                self.alignment_to_sequence.append(-1)

            else:

                self.alignment_to_sequence.append(sequence_position)

                self.sequence_to_alignment.append(alignment_position)

                sequence_position += 1

    @property
    def alignment_length(self) -> int:

        return len(self.aligned)

    @property
    def sequence_length(self) -> int:

        return len(self.ungapped)

    def alignment_to_seq(self, position: int) -> int:

        if position < 0 or position >= len(self.alignment_to_sequence):
            raise IndexError("Alignment position out of range.")

        return self.alignment_to_sequence[position]

    def seq_to_alignment(self, position: int) -> int:

        if position < 0 or position >= len(self.sequence_to_alignment):
            raise IndexError("Sequence position out of range.")

        return self.sequence_to_alignment[position]

    def base_at_alignment(self, position: int) -> str:

        return self.aligned[position]

    def base_at_sequence(self, position: int) -> str:

        return self.ungapped[position]
