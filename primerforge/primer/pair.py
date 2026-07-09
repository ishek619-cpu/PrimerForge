"""
Primer pair generation.
"""

from primerforge.models.pair import PrimerPair
from primerforge.models.primer import Primer


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

                if reverse.start <= forward.end:
                    continue

                product = reverse.end - forward.start + 1

                if product < self.min_product:
                    continue

                if product > self.max_product:
                    continue

                pairs.append(
                    PrimerPair(
                        forward=forward,
                        reverse=reverse,
                        product_size=product,
                    )
                )

        return pairs
