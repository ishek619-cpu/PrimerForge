"""
Candidate search.

Uses a seed index to rapidly identify candidate
background sequences for primer matching.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class CandidateHit:

    sequence_id: int

    position: int

    strand: str

    seed: str


class CandidateSearcher:

    def __init__(

        self,

        seed_length: int = 18,

    ):

        self.seed_length = seed_length

    ############################################################

    def reverse_complement(

        self,

        sequence: str,

    ) -> str:

        table = str.maketrans(

            "ACGTN",

            "TGCAN",

        )

        return sequence.translate(table)[::-1]

    ############################################################

    def search(

        self,

        primer: str,

        index,

    ) -> list[CandidateHit]:

        primer = primer.upper()

        seed = primer[-self.seed_length:]

        hits = []

        ########################################################
        # Forward strand
        ########################################################

        for sequence_id, position in index.lookup(seed):

            hits.append(

                CandidateHit(

                    sequence_id=sequence_id,

                    position=position,

                    strand="+",

                    seed=seed,

                )

            )

        ########################################################
        # Reverse strand
        ########################################################

        rc = self.reverse_complement(seed)

        for sequence_id, position in index.lookup(rc):

            hits.append(

                CandidateHit(

                    sequence_id=sequence_id,

                    position=position,

                    strand="-",

                    seed=rc,

                )

            )

        return hits
