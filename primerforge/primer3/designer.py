"""
Primer3 interface.
"""

from pathlib import Path
from typing import List

import primer3

from primerforge.models.primer import Primer


class Primer3Designer:

    def design(
        self,
        template: str,
        product_size=(80, 150),
    ) -> List[Primer]:

        result = primer3.bindings.designPrimers(
            {
                "SEQUENCE_TEMPLATE": template,
            },
            {
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 25,

                "PRIMER_OPT_TM": 60.0,
                "PRIMER_MIN_TM": 58.0,
                "PRIMER_MAX_TM": 62.0,

                "PRIMER_MIN_GC": 40,
                "PRIMER_MAX_GC": 60,

                "PRIMER_PRODUCT_SIZE_RANGE":
                    [list(product_size)],

                "PRIMER_NUM_RETURN": 20,
            },
        )

        primers = []

        count = result["PRIMER_LEFT_NUM_RETURNED"]

        for i in range(count):

            seq = result[f"PRIMER_LEFT_{i}_SEQUENCE"]

            start, length = result[f"PRIMER_LEFT_{i}"]

            tm = result[f"PRIMER_LEFT_{i}_TM"]

            gc = result[f"PRIMER_LEFT_{i}_GC_PERCENT"]

            primers.append(

                Primer(
                    sequence=seq,
                    start=start + 1,
                    end=start + length,
                    strand="+",
                    length=length,
                    tm=tm,
                    gc=gc,
                )

            )

        return primers
