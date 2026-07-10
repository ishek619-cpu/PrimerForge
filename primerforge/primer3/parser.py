"""
Primer3 result parser.
"""

from primerforge.models.primer import Primer
from primerforge.models.pair import PrimerPair


class Primer3Parser:

    def parse_left(self, result):

        primers = []

        n = result.get(
            "PRIMER_LEFT_NUM_RETURNED",
            0,
        )

        for i in range(n):

            start, length = result[
                f"PRIMER_LEFT_{i}"
            ]

            primers.append(

                Primer(
                    sequence=result[
                        f"PRIMER_LEFT_{i}_SEQUENCE"
                    ],
                    start=start + 1,
                    end=start + length,
                    strand="+",
                    length=length,
                    tm=result[
                        f"PRIMER_LEFT_{i}_TM"
                    ],
                    gc=result[
                        f"PRIMER_LEFT_{i}_GC_PERCENT"
                    ],
                )

            )

        return primers

    def parse_right(self, result):

        primers = []

        n = result.get(
            "PRIMER_RIGHT_NUM_RETURNED",
            0,
        )

        for i in range(n):

            start, length = result[
                f"PRIMER_RIGHT_{i}"
            ]

            primers.append(

                Primer(
                    sequence=result[
                        f"PRIMER_RIGHT_{i}_SEQUENCE"
                    ],
                    start=start + 1,
                    end=start + length,
                    strand="-",
                    length=length,
                    tm=result[
                        f"PRIMER_RIGHT_{i}_TM"
                    ],
                    gc=result[
                        f"PRIMER_RIGHT_{i}_GC_PERCENT"
                    ],
                )

            )

        return primers

    def parse_pairs(self, result):

        forward = self.parse_left(result)

        reverse = self.parse_right(result)

        pairs = []

        n = min(
            len(forward),
            len(reverse),
        )

        for i in range(n):

            product = result.get(
                f"PRIMER_PAIR_{i}_PRODUCT_SIZE",
                reverse[i].end - forward[i].start + 1,
            )

            pair = PrimerPair(
                forward=forward[i],
                reverse=reverse[i],
                product_size=product,
            )

            if (
                f"PRIMER_PAIR_{i}_PENALTY"
                in result
            ):
                pair.breakdown[
                    "primer3_penalty"
                ] = round(
                    result[
                        f"PRIMER_PAIR_{i}_PENALTY"
                    ],
                    3,
                )

            pairs.append(pair)

        return pairs
