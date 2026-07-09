"""
JSON report generator.
"""

import json
from pathlib import Path

from primerforge.models.pair import PrimerPair


class JSONReport:
    """
    Generate JSON report.
    """

    def write(
        self,
        pairs: list[PrimerPair],
        output: Path,
    ):

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        report = []

        for pair in pairs:

            report.append({

                "forward": pair.forward.sequence,

                "reverse": pair.reverse.sequence,

                "forward_tm": round(
                    pair.forward.tm,
                    2,
                ),

                "reverse_tm": round(
                    pair.reverse.tm,
                    2,
                ),

                "forward_gc": pair.forward.gc,

                "reverse_gc": pair.reverse.gc,

                "product_size": pair.product_size,

                "score": pair.score,

            })

        output.write_text(

            json.dumps(
                report,
                indent=4,
            ),

            encoding="utf-8",

        )

