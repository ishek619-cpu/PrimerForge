"""
Predict PCR products from paired BLAST hits.
"""

from primerforge.models.pcr_product import PCRProduct
from primerforge.specificity.models import BlastHit


class PCRProductPredictor:
    """
    Predict PCR amplification products from a pair of BLAST hits.
    """

    def __init__(
        self,
        min_product: int = 50,
        max_product: int = 1000,
    ):

        self.min_product = min_product

        self.max_product = max_product

    def predict(
        self,
        forward_hit: BlastHit,
        reverse_hit: BlastHit,
    ) -> PCRProduct | None:

        #
        # Primers must align to the same reference.
        #
        if forward_hit.accession != reverse_hit.accession:

            return None

        #
        # Forward primer must be upstream.
        #
        if forward_hit.sstart >= reverse_hit.sstart:

            return None

        #
        # Calculate amplicon size.
        #
        product_size = (

            reverse_hit.send

            - forward_hit.sstart

            + 1

        )

        passed = (

            self.min_product

            <= product_size

            <= self.max_product

        )

        reason = None

        if product_size < self.min_product:

            reason = "Product too small"

        elif product_size > self.max_product:

            reason = "Product too large"

        return PCRProduct(

            chromosome=forward_hit.accession,

            forward_start=forward_hit.sstart,

            reverse_end=reverse_hit.send,

            product_size=product_size,

            forward_strand=forward_hit.strand,

            reverse_strand=reverse_hit.strand,

            passed=passed,

            reason=reason,

        )
