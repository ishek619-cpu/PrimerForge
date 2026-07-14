"""
Alignment model for EnvoPrimer.

Owns all aligned sequences and provides alignment-wide utilities.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio import AlignIO

from envoprimer_v2.models.sequence import Sequence


@dataclass(slots=True)
class Alignment:

    sequences: list[Sequence] = field(default_factory=list)

    @classmethod
    def from_fasta(cls, fasta):

        msa = AlignIO.read(
            fasta,
            "fasta",
        )

        sequences = []

        for record in msa:

            header = record.description

            accession = record.id

            species = "Unknown"

            words = header.split()

            if len(words) >= 2:

                species = words[1]

                if len(words) >= 3:

                    species += " " + words[2]

            sequences.append(

                Sequence(

                    id=record.id,

                    accession=accession,

                    species=species,

                    aligned=str(record.seq),

                )

            )

        return cls(
            sequences=sequences,
        )

    @property
    def length(self):

        if not self.sequences:

            return 0

        return self.sequences[0].alignment_length

    @property
    def count(self):

        return len(
            self.sequences,
        )

    def column(self, position):

        return [

            sequence.base_at_alignment(position)

            for sequence in self.sequences

        ]

    def consensus(self):

        consensus = []

        for column in range(self.length):

            bases = [

                b

                for b in self.column(column)

                if b != "-"

            ]

            if not bases:

                consensus.append("-")

                continue

            consensus.append(

                max(

                    set(bases),

                    key=bases.count,

                )

            )

        return "".join(consensus)

    def get(self, accession):

        for sequence in self.sequences:

            if sequence.accession == accession:

                return sequence

        raise KeyError(accession)
