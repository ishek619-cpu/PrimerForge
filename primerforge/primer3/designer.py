"""
Primer3 design engine.
"""

from __future__ import annotations

import primer3

from primerforge.models.region import Region
from primerforge.primer3.parser import Primer3Parser
from primerforge.primer3.settings import PRIMER3_SETTINGS
from primerforge.reference.annotation import Gene


class Primer3Designer:
    """
    Interface to Primer3.
    """

    def __init__(self):

        self.parser = Primer3Parser()

    def design(
        self,
        template: str,
        product_size: tuple[int, int] = (80, 150),
        region: Region | Gene | None = None,
    ):
        """
        Design primers for a DNA template.

        Parameters
        ----------
        template
            DNA sequence.

        product_size
            Desired PCR product size range.

        region
            Either a Region object (conservation analysis)
            or a Gene object (GenBank annotation).
        """

        settings = PRIMER3_SETTINGS.copy()

        settings[
            "PRIMER_PRODUCT_SIZE_RANGE"
        ] = [
            list(product_size)
        ]

        seq_args = {

            "SEQUENCE_TEMPLATE": template,

        }

        #
        # Restrict Primer3 to a selected region.
        #
        if region is not None:

            #
            # Gene annotation
            #
            if isinstance(
                region,
                Gene,
            ):

                start = region.start - 1

                length = region.length

            #
            # Conserved region
            #
            else:

                start = region.start

                length = region.length

            seq_args[
                "SEQUENCE_INCLUDED_REGION"
            ] = [
                start,
                length,
            ]

        result = primer3.bindings.design_primers(

            seq_args=seq_args,

            global_args=settings,

        )

        forward = self.parser.parse_left(
            result,
        )

        reverse = self.parser.parse_right(
            result,
        )

        return (
            forward,
            reverse,
        )

    def design_gene(
        self,
        template: str,
        gene: Gene,
        product_size: tuple[int, int] = (80, 150),
    ):
        """
        Convenience wrapper for designing primers
        directly from a Gene annotation.
        """

        return self.design(

            template=template,

            product_size=product_size,

            region=gene,

        )
