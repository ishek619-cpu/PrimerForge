"""
Primer3 result parser.

Converts Primer3 output into EnvoPrimer Primer objects while
preserving absolute reference coordinates.
"""

from __future__ import annotations

from envoprimer_v2.models.primer import Primer


class PrimerParser:

    def parse(self, result: dict, window_start: int):

        forward = []
        reverse = []

        n_left = result.get(
            "PRIMER_LEFT_NUM_RETURNED",
            0,
        )

        n_right = result.get(
            "PRIMER_RIGHT_NUM_RETURNED",
            0,
        )

        #
        # Forward primers
        #
        for i in range(n_left):

            start, length = result[
                f"PRIMER_LEFT_{i}"
            ]

            absolute_start = window_start + start + 1

            absolute_end = (
                absolute_start
                + length
                - 1
            )

            penalty = result.get(
                f"PRIMER_LEFT_{i}_PENALTY",
                0.0,
            )

            forward.append(

                Primer(

                    sequence=result[
                        f"PRIMER_LEFT_{i}_SEQUENCE"
                    ],

                    start=absolute_start,

                    end=absolute_end,

                    strand="+",

                    tm=result[
                        f"PRIMER_LEFT_{i}_TM"
                    ],

                    gc=result[
                        f"PRIMER_LEFT_{i}_GC_PERCENT"
                    ],

                    length=length,

                    penalty=penalty,

                )

            )

        #
        # Reverse primers
        #
        for i in range(n_right):

            start, length = result[
                f"PRIMER_RIGHT_{i}"
            ]

            absolute_start = window_start + start + 1

            absolute_end = (
                absolute_start
                + length
                - 1
            )

            penalty = result.get(
                f"PRIMER_RIGHT_{i}_PENALTY",
                0.0,
            )

            reverse.append(

                Primer(

                    sequence=result[
                        f"PRIMER_RIGHT_{i}_SEQUENCE"
                    ],

                    start=absolute_start,

                    end=absolute_end,

                    strand="-",

                    tm=result[
                        f"PRIMER_RIGHT_{i}_TM"
                    ],

                    gc=result[
                        f"PRIMER_RIGHT_{i}_GC_PERCENT"
                    ],

                    length=length,

                    penalty=penalty,

                )

            )

        return forward, reverse
