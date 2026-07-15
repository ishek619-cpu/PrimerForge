"""
Binding location reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class BindingLocation:

    sequence_id: int

    start: int

    end: int

    strand: str


class BindingLocator:

    def __init__(

        self,

        seed_length: int = 18,

    ):

        self.seed_length = seed_length

    ############################################################

    def locate(

        self,

        primer: str,

        candidate,

    ) -> BindingLocation:

        primer_length = len(primer)

        offset = (

            primer_length

            - self.seed_length

        )

        ########################################################
        # Forward strand
        ########################################################

        if candidate.strand == "+":

            start = candidate.position - offset

            end = start + primer_length

        ########################################################
        # Reverse strand
        ########################################################

        else:

            #
            # Candidate position is the FIRST base
            # of the reverse-complement seed.
            #
            # Convert to the FULL primer coordinates.
            #

            start = candidate.position - offset

            end = start + primer_length

        return BindingLocation(

            sequence_id=candidate.sequence_id,

            start=start,

            end=end,

            strand=candidate.strand,

        )
