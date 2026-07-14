"""
Primer pair generation.
"""

from envoprimer.models.pair import PrimerPair
from envoprimer.models.primer import Primer


class PrimerPairGenerator:
    """
    Generate primer pairs from candidate primers.
    """

    def __init__(
        self,
        min_product: int = 80,
        max_product: int = 150,
    ):

        self.min_product = min_product
        self.max_product = max_product

    def generate(
        self,
        forward_primers: list[Primer],
        reverse_primers: list[Primer],
    ) -> list[PrimerPair]:

        pairs = []

        for forward in forward_primers:

            for reverse in reverse_primers:

                #
                # Use universal coordinates
                #
                forward_start = (
                    forward.coordinate.start
                )

                forward_end = (
                    forward.coordinate.end
                )

                reverse_start = (
                    reverse.coordinate.start
                )

                reverse_end = (
                    reverse.coordinate.end
                )

                #
                # Correct orientation
                #
                if reverse_start <= forward_end:

                    continue

                #
                # Product size
                #
                product = (

                    reverse_end

                    - forward_start

                    + 1

                )

                #
                # Product size filter
                #
                if product < self.min_product:

                    continue

                if product > self.max_product:

                    continue

                #
                # Primer Tm compatibility
                #
                if abs(

                    forward.tm

                    - reverse.tm

                ) > 2.0:

                    continue

                #
                # Primer GC compatibility
                #
                if abs(

                    forward.gc

                    - reverse.gc

                ) > 15.0:

                    continue

                pairs.append(

                    PrimerPair(

                        forward=forward,

                        reverse=reverse,

                        product_size=product,

                    )

                )

        #
        # Best-balanced primer pairs first
        #
        pairs.sort(

            key=lambda pair: (

                abs(
                    pair.forward.tm
                    - pair.reverse.tm
                ),

                abs(
                    pair.forward.gc
                    - pair.reverse.gc
                ),

                abs(
                    pair.product_size
                    - 120
                ),

            )

        )

        return pairs
