"""
PCR product prediction.
"""

from __future__ import annotations

from envoprimer_v2.pcr.models import BindingSite, PCRProduct


class ProductPredictor:
    """
    Predict PCR amplicons from accepted forward and reverse bindings.
    """

    def __init__(
        self,
        min_size: int = 50,
        max_size: int = 2000,
    ):
        self.min_size = min_size
        self.max_size = max_size

    def predict(
        self,
        forward_sites: list[BindingSite],
        reverse_sites: list[BindingSite],
    ) -> list[PCRProduct]:

        products = []

        for f in forward_sites:

            for r in reverse_sites:

                if f.sequence_id != r.sequence_id:
                    continue

                if f.start >= r.start:
                    continue

                size = r.end - f.start

                if size < self.min_size:
                    continue

                if size > self.max_size:
                    continue

                products.append(
                    PCRProduct(
                        sequence_id=f.sequence_id,
                        forward_site=f,
                        reverse_site=r,
                        start=f.start,
                        end=r.end,
                        size=size,
                        predicted=True,
                    )
                )

        return products
