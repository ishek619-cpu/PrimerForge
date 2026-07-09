"""
Primer3 result parser.
"""

from primerforge.models.primer import Primer


class Primer3Parser:

    def parse_left(self, result):

        primers = []

        n = result.get("PRIMER_LEFT_NUM_RETURNED", 0)

        for i in range(n):

            start, length = result[f"PRIMER_LEFT_{i}"]

            primers.append(

                Primer(
                    sequence=result[f"PRIMER_LEFT_{i}_SEQUENCE"],
                    start=start + 1,
                    end=start + length,
                    strand="+",
                    length=length,
                    tm=result[f"PRIMER_LEFT_{i}_TM"],
                    gc=result[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                )

            )

        return primers

    def parse_right(self, result):

        primers = []

        n = result.get("PRIMER_RIGHT_NUM_RETURNED", 0)

        for i in range(n):

            start, length = result[f"PRIMER_RIGHT_{i}"]

            primers.append(

                Primer(
                    sequence=result[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                    start=start + 1,
                    end=start + length,
                    strand="-",
                    length=length,
                    tm=result[f"PRIMER_RIGHT_{i}_TM"],
                    gc=result[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                )

            )

        return primers
