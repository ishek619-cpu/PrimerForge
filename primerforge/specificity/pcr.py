"""
PCR product prediction.
"""

from dataclasses import dataclass


@dataclass(slots=True)
class PCRProduct:

    chromosome: str

    forward_start: int

    forward_end: int

    reverse_start: int

    reverse_end: int

    strand: str

    size: int

    identity: float


class PCRProductFinder:
    """
    Predict PCR products from BLAST hits.
    """

    def find_products(
        self,
        forward_hits,
        reverse_hits,
        min_size: int = 50,
        max_size: int = 1000,
    ):

        products = []

        for f in forward_hits:

            for r in reverse_hits:

                if f["subject"] != r["subject"]:
                    continue

                if f["strand"] != "plus":
                    continue

                if r["strand"] != "minus":
                    continue

                start = min(
                    f["subject_start"],
                    f["subject_end"],
                )

                end = max(
                    r["subject_start"],
                    r["subject_end"],
                )

                if end <= start:
                    continue

                size = end - start + 1

                if size < min_size:
                    continue

                if size > max_size:
                    continue

                products.append(

                    PCRProduct(

                        chromosome=f["subject"],

                        forward_start=f["subject_start"],

                        forward_end=f["subject_end"],

                        reverse_start=r["subject_start"],

                        reverse_end=r["subject_end"],

                        strand="plus",

                        size=size,

                        identity=min(
                            f["identity"],
                            r["identity"],
                        ),

                    )

                )

        products.sort(
            key=lambda x: x.size,
        )

        return products
