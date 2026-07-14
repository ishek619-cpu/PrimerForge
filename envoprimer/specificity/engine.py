"""
High-performance species specificity engine.

Uses:
BackgroundIndex
        ↓
CandidateSearcher
        ↓
PrimerAlignment
        ↓
PCREvaluator

instead of scanning every sequence for every primer.
"""

from __future__ import annotations

from envoprimer.specificity.index import BackgroundIndex
from envoprimer.specificity.pcr import PCREvaluator


class SpecificityEngine:

    def __init__(

        self,

        seed_length: int = 18,

    ):

        self.index = BackgroundIndex(

            k=seed_length,

        )

        self.pcr = PCREvaluator()

        self._built = False

    ############################################################

    def build(

        self,

        background_sequences: list[str],

    ):

        if self._built:

            return

        print()

        print("=" * 60)
        print("BUILDING BACKGROUND INDEX")
        print("=" * 60)

        self.index.build(

            background_sequences,

        )

        self._built = True

    ############################################################

    def score(

        self,

        pair,

        background_sequences,

    ):

        ########################################################
        # Build index only once
        ########################################################

        self.build(

            background_sequences,

        )

        ########################################################
        # Evaluate possible PCR products
        ########################################################

        hits = self.pcr.evaluate(

            pair,

            self.index,

        )

        ########################################################

        n = len(hits)

        ########################################################
        # Weighted penalty
        ########################################################

        if n == 0:

            score = 100.0

        elif n == 1:

            score = 85.0

        elif n <= 3:

            score = 60.0

        elif n <= 10:

            score = 30.0

        else:

            score = 0.0

        pair.specificity_score = score

        pair.metadata["background_hits"] = n

        pair.metadata["background_details"] = [

            {

                "sequence_id": hit.sequence_id,

                "product_size": hit.product_size,

                "forward_score": hit.forward_score,

                "reverse_score": hit.reverse_score,

            }

            for hit in hits

        ]

        return score
