"""
Primer3 design engine.
"""

import primer3

from primerforge.models.region import Region
from primerforge.primer3.parser import Primer3Parser
from primerforge.primer3.settings import PRIMER3_SETTINGS


class Primer3Designer:
    """
    Interface to Primer3.
    """

    def __init__(self):

        self.parser = Primer3Parser()

    def design(
        self,
        template: str,
        product_size=(80, 150),
        region: Region | None = None,
    ):
        """
        Design primers for a DNA template.

        If a conserved region is supplied, Primer3 searches
        only within that region. Otherwise the whole template
        is searched.
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

        if region is not None:

            seq_args[
                "SEQUENCE_INCLUDED_REGION"
            ] = [
                region.start,
                region.length,
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
