"""
PCR product prediction.
"""

from dataclasses import dataclass


@dataclass(slots=True)
class PCRProduct:

    forward_start: int

    reverse_end: int

    size: int

    chromosome: str

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

                if (
                    f["subject"] !=
                    r["subject"]
                ):
                    continue

                if (
                    f["query_start"] >=
                    r["query_start"]
                ):
                    continue

                size = (

                    r["query_end"]

                    -

                    f["query_start"]

                    +

                    1

                )

                if (
                    size < min_size
                    or
                    size > max_size
                ):
                    continue

                products.append(

                    PCRProduct(

                        forward_start=f["query_start"],

                        reverse_end=r["query_end"],

                        size=size,

                        chromosome=f["subject"],

                        identity=min(

                            f["identity"],

                            r["identity"],

                        ),

                    )

                )

        return products
