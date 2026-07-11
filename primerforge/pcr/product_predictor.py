"""
Predict PCR products from paired BLAST hits.
"""

from primerforge.models.pcr_product import PCRProduct
from primerforge.reference.coordinates import Coordinate
from primerforge.specificity.models import BlastHit


class PCRProductPredictor:
    """
    Predict PCR amplification products from paired BLAST hits.
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
        # Both primers must hit the same sequence.
        #
        if forward_hit.accession != reverse_hit.accession:

            return None

        #
        # Forward primer must occur before reverse primer.
        #
        if forward_hit.sstart >= reverse_hit.sstart:

            return None

        #
        # Represent primer positions using Coordinates.
        #
        forward = Coordinate(
            start=forward_hit.sstart,
            end=forward_hit.send,
            system="reference",
        )

        reverse = Coordinate(
            start=reverse_hit.sstart,
            end=reverse_hit.send,
            system="reference",
        )

        #
        # Calculate amplicon size.
        #
        product_size = (
            reverse.end
            - forward.start
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

            forward_start=forward.start,

            reverse_end=reverse.end,

            product_size=product_size,

            forward_strand=forward_hit.strand,

            reverse_strand=reverse_hit.strand,

            passed=passed,

            reason=reason,

        )
