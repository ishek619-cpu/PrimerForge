"""
Primer pair generation.

EnvoPrimer v2
"""

from __future__ import annotations

from envoprimer_v2.models.pair import PrimerPair
from envoprimer_v2.models.primer import Primer


class PrimerPairGenerator:

    def __init__(

        self,

        minimum_product: int = 70,

        maximum_product: int = 200,

        maximum_tm_difference: float = 2.0,

        maximum_gc_difference: float = 15.0,

    ):

        self.minimum_product = minimum_product

        self.maximum_product = maximum_product

        self.maximum_tm_difference = maximum_tm_difference

        self.maximum_gc_difference = maximum_gc_difference

    def generate(

        self,

        forward_primers: list[Primer],

        reverse_primers: list[Primer],

    ) -> list[PrimerPair]:

        pairs = []

        for forward in forward_primers:

            for reverse in reverse_primers:

                #
                # Orientation
                #

                if reverse.start <= forward.end:

                    continue

                #
                # Product size
                #

                product = (

                    reverse.end

                    - forward.start

                    + 1

                )

                if product < self.minimum_product:

                    continue

                if product > self.maximum_product:

                    continue

                #
                # Tm compatibility
                #

                tm_difference = abs(

                    forward.tm

                    - reverse.tm

                )

                if tm_difference > self.maximum_tm_difference:

                    continue

                #
                # GC compatibility
                #

                gc_difference = abs(

                    forward.gc

                    - reverse.gc

                )

                if gc_difference > self.maximum_gc_difference:

                    continue

                penalty = (

                    forward.penalty

                    + reverse.penalty

                )

                pair = PrimerPair(

                    forward=forward,

                    reverse=reverse,

                    product_size=product,

                    primer3_penalty=round(

                        penalty,

                        3,

                    ),

                )

                pairs.append(

                    pair,

                )

        pairs.sort(

            key=lambda pair: (

                pair.primer3_penalty,

                abs(pair.product_size - 130),

                pair.tm_difference,

                pair.gc_difference,

            )

        )

        return pairs
