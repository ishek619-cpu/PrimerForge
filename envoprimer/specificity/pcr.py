"""
PCR specificity evaluator.

Combines candidate search and primer alignment to estimate
whether a primer pair can amplify a background sequence.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer.specificity.search import CandidateSearcher
from envoprimer.specificity.alignment import PrimerAlignment


@dataclass(slots=True)
class PCRHit:

    sequence_id: int

    forward_score: float

    reverse_score: float

    product_size: int

    passed: bool


class PCREvaluator:

    def __init__(

        self,

        minimum_product: int = 70,

        maximum_product: int = 200,

    ):

        self.search = CandidateSearcher()

        self.align = PrimerAlignment()

        self.minimum_product = minimum_product

        self.maximum_product = maximum_product

    ####################################################################

    def evaluate(

        self,

        pair,

        index,

    ) -> list[PCRHit]:

        hits = []

        forward_hits = self.search.search(

            pair.forward.sequence,

            index,

        )

        reverse_hits = self.search.search(

            pair.reverse.sequence,

            index,

        )

        ############################################################

        reverse_lookup = {}

        for hit in reverse_hits:

            reverse_lookup.setdefault(

                hit.sequence_id,

                [],

            ).append(hit)

        ############################################################

        for forward in forward_hits:

            if forward.sequence_id not in reverse_lookup:

                continue

            sequence = index.sequence(

                forward.sequence_id,

            )

            f = self.align.align(

                pair.forward.sequence,

                sequence,

                forward.position,

                forward.strand,

            )

            if not f.passed:

                continue

            for reverse in reverse_lookup[

                forward.sequence_id

            ]:

                r = self.align.align(

                    pair.reverse.sequence,

                    sequence,

                    reverse.position,

                    reverse.strand,

                )

                if not r.passed:

                    continue

                ################################################

                product = (

                    reverse.position

                    - forward.position

                    + len(pair.reverse.sequence)

                )

                if product < self.minimum_product:

                    continue

                if product > self.maximum_product:

                    continue

                hits.append(

                    PCRHit(

                        sequence_id=forward.sequence_id,

                        forward_score=f.score,

                        reverse_score=r.score,

                        product_size=product,

                        passed=True,

                    )

                )

        return hits
