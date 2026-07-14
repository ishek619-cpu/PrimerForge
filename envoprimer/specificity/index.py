"""
Background sequence index.

Builds a positional k-mer index for rapid primer
specificity testing.
"""

from __future__ import annotations

from collections import defaultdict


class BackgroundIndex:

    def __init__(

        self,

        k: int = 18,

    ):

        self.k = k

        #
        # kmer -> [(sequence_id, position)]
        #
        self.index = defaultdict(list)

        self.sequences = []

    ############################################################

    def build(

        self,

        sequences: list[str],

    ):

        self.index.clear()

        self.sequences = [

            seq.upper()

            for seq in sequences

        ]

        total = 0

        for sequence_id, sequence in enumerate(

            self.sequences

        ):

            if len(sequence) < self.k:

                continue

            for position in range(

                len(sequence) - self.k + 1

            ):

                kmer = sequence[

                    position:position + self.k

                ]

                if "-" in kmer:

                    continue

                if "N" in kmer:

                    continue

                self.index[kmer].append(

                    (

                        sequence_id,

                        position,

                    )

                )

                total += 1

        print()

        print("=" * 60)
        print("BACKGROUND INDEX")
        print("=" * 60)

        print(f"Sequences      : {len(self.sequences):,}")
        print(f"K-mer length   : {self.k}")
        print(f"Unique k-mers  : {len(self.index):,}")
        print(f"Indexed k-mers : {total:,}")

    ############################################################

    def lookup(

        self,

        kmer: str,

    ):

        return self.index.get(

            kmer.upper(),

            [],

        )

    ############################################################

    def sequence(

        self,

        sequence_id: int,

    ) -> str:

        return self.sequences[sequence_id]
