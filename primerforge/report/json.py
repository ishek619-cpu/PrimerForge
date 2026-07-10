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

        for rank, pair in enumerate(
            pairs,
            start=1,
        ):

            report.append(

                {

                    "rank": rank,

                    "score": pair.score,

                    "forward": {

                        "sequence": pair.forward.sequence,

                        "tm": round(
                            pair.forward.tm,
                            2,
                        ),

                        "gc": round(
                            pair.forward.gc,
                            2,
                        ),

                    },

                    "reverse": {

                        "sequence": pair.reverse.sequence,

                        "tm": round(
                            pair.reverse.tm,
                            2,
                        ),

                        "gc": round(
                            pair.reverse.gc,
                            2,
                        ),

                    },

                    "product_size": pair.product_size,

                    "population": {

                        "conservation": round(
                            pair.population_conservation,
                            2,
                        ),

                        "coverage": round(
                            pair.population_coverage,
                            2,
                        ),

                        "details": pair.population_result,

                    },

                    "specificity": {

                        "score": round(
                            pair.specificity_score,
                            2,
                        ),

                        "passed": pair.passed_specificity,

                        "details": (
                            None
                            if pair.specificity_result is None
                            else {

                                "target_hits": len(
                                    pair.specificity_result.target_hits
                                ),

                                "off_target_hits": len(
                                    pair.specificity_result.off_target_hits
                                ),

                                "rejection_reason": (
                                    pair.specificity_result.rejection_reason
                                ),

                            }
                        ),

                    },

                    "score_breakdown": pair.breakdown,

                }

            )

        output.write_text(

            json.dumps(
                report,
                indent=4,
            ),

            encoding="utf-8",

        )

        return output
