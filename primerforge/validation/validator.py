"""
Complete primer validation engine.
"""

from pathlib import Path

from primerforge.models.pair import PrimerPair
from primerforge.models.snp import SNP

from primerforge.primer3.thermo import ThermoAnalyzer
from primerforge.validation.conservation import PrimerConservation
from primerforge.primer.snp_optimizer import SNPOptimizer


class PrimerValidator:
    """
    Validate primer pairs using thermodynamics,
    conservation and SNP analysis.
    """

    def __init__(self):

        self.thermo = ThermoAnalyzer()
        self.conservation = PrimerConservation()
        self.snp = SNPOptimizer()

    def validate(
        self,
        pair: PrimerPair,
        alignment: Path,
        snps: list[SNP],
    ):

        pair.forward = self.thermo.evaluate(pair.forward)
        pair.reverse = self.thermo.evaluate(pair.reverse)

        heterodimer = self.thermo.heterodimer(
            pair.forward,
            pair.reverse,
        )

        forward = self.conservation.evaluate(
            pair.forward,
            alignment,
        )

        reverse = self.conservation.evaluate(
            pair.reverse,
            alignment,
        )

        forward_snp = self.snp.score(
            pair.forward,
            snps,
        )

        reverse_snp = self.snp.score(
            pair.reverse,
            snps,
        )

        thermo = 100.0

        thermo -= pair.forward.hairpin_score
        thermo -= pair.reverse.hairpin_score

        thermo -= pair.forward.self_dimer_score
        thermo -= pair.reverse.self_dimer_score

        thermo -= heterodimer

        thermo = max(
            thermo,
            0.0,
        )

        conservation = (
            forward["coverage"] +
            reverse["coverage"]
        ) / 2.0

        snp = (
            forward_snp +
            reverse_snp
        ) / 2.0

        final = (
            thermo * 0.40 +
            conservation * 0.40 +
            snp * 0.20
        )

        if final >= 90:
            status = "PASS"
        elif final >= 75:
            status = "WARNING"
        else:
            status = "FAIL"

        return {

            "thermodynamics": round(
                thermo,
                2,
            ),

            "forward_conservation": forward,

            "reverse_conservation": reverse,

            "snp_score": round(
                snp,
                2,
            ),

            "final_score": round(
                final,
                2,
            ),

            "status": status,

        }
