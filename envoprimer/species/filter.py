"""
Sequence quality filtering.
"""

from __future__ import annotations

from collections import OrderedDict

from Bio.SeqRecord import SeqRecord


class SequenceFilter:

    def __init__(

        self,

        minimum_length: int = 500,

    ):

        self.minimum_length = minimum_length

    def clean(

        self,

        records: list[SeqRecord],

    ) -> list[SeqRecord]:

        unique = OrderedDict()

        for record in records:

            sequence = str(record.seq).upper()

            #
            # Remove Ns
            #
            if sequence.count("N") > 5:

                continue

            #
            # Remove short sequences
            #
            if len(sequence) < self.minimum_length:

                continue

            #
            # Remove duplicates
            #
            if sequence in unique:

                continue

            unique[sequence] = record

        return list(unique.values())
