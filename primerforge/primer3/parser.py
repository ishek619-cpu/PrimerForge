"""
Primer3 result parser.
"""

from primerforge.models.pair import PrimerPair
from primerforge.models.primer import Primer
from primerforge.reference.coordinates import Coordinate


class Primer3Parser:
    """
    Parse Primer3 output into PrimerForge models.
    """

    def parse_left(
        self,
        result,
    ):

        primers = []

        n = result.get(
            "PRIMER_LEFT_NUM_RETURNED",
            0,
        )

        for i in range(n):

            start, length = result[
                f"PRIMER_LEFT_{i}"
            ]

            start += 1

            end = start + length - 1

            primers.append(

                Primer(

                    sequence=result[
                        f"PRIMER_LEFT_{i}_SEQUENCE"
                    ],

                    start=start,

                    end=end,

                    strand="+",

                    length=length,

                    tm=result[
                        f"PRIMER_LEFT_{i}_TM"
                    ],

                    gc=result[
                        f"PRIMER_LEFT_{i}_GC_PERCENT"
                    ],

                    coordinate=Coordinate(

                        start=start,

                        end=end,

                        system="reference",

                    ),

                )

            )

        return primers

    def parse_right(
        self,
        result,
    ):

        primers = []

        n = result.get(
            "PRIMER_RIGHT_NUM_RETURNED",
            0,
        )

        for i in range(n):

            start, length = result[
                f"PRIMER_RIGHT_{i}"
            ]

            start += 1

            end = start + length - 1

            primers.append(

                Primer(

                    sequence=result[
                        f"PRIMER_RIGHT_{i}_SEQUENCE"
                    ],

                    start=start,

                    end=end,

                    strand="-",

                    length=length,

                    tm=result[
                        f"PRIMER_RIGHT_{i}_TM"
                    ],

                    gc=result[
                        f"PRIMER_RIGHT_{i}_GC_PERCENT"
                    ],

                    coordinate=Coordinate(

                        start=start,

                        end=end,

                        system="reference",

                    ),

                )

            )

        return primers

    def parse_pairs(
        self,
        result,
    ):

        forward = self.parse_left(
            result,
        )

        reverse = self.parse_right(
            result,
        )

        pairs = []

        n = min(
            len(forward),
            len(reverse),
        )

        for i in range(n):

            product = result.get(

                f"PRIMER_PAIR_{i}_PRODUCT_SIZE",

                reverse[i].coordinate.end
                - forward[i].coordinate.start
                + 1,

            )

            pair = PrimerPair(

                forward=forward[i],

                reverse=reverse[i],

                product_size=product,

            )

            penalty = result.get(
                f"PRIMER_PAIR_{i}_PENALTY",
            )

            if penalty is not None:

                pair.breakdown[
                    "primer3_penalty"
                ] = round(
                    penalty,
                    3,
                )

            pairs.append(
                pair,
            )

        return pairs
