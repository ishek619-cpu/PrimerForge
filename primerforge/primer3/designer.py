"""
Primer3 design engine.
"""

import primer3

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
    ):
        """
        Design primers for a DNA template.
        """

        settings = PRIMER3_SETTINGS.copy()

        settings["PRIMER_PRODUCT_SIZE_RANGE"] = [
            list(product_size)
        ]

        result = primer3.bindings.design_primers(
            seq_args={
                "SEQUENCE_TEMPLATE": template,
            },
            global_args=settings,
        )

        forward = self.parser.parse_left(result)

        reverse = self.parser.parse_right(result)

        return forward, reverse
