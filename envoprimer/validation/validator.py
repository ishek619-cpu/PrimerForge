"""
Primer validation.
"""

from pathlib import Path

from envoprimer.analysis.snps import SNPFinder
from envoprimer.models.pair import PrimerPair
from envoprimer.models.snp import SNP
from envoprimer.primer3.thermo import ThermoAnalyzer
from envoprimer.validation.conservation import PrimerConservation


class SNPScorer:
    """
    Score primers against known SNPs.
    """

    def score(
        self,
        primer,
        snps: list[SNP],
    ) -> float:

        if not snps:
            return 100.0

        overlapping = 0

        for snp in snps:

            if primer.start <= snp.position <= primer.end:
                overlapping += 1

        return round(
            max(
                0.0,
                100.0 - overlapping * 10.0,
            ),
            2,
        )


class PrimerValidator:
    """
    Validate primer pairs.
    """

    def __init__(self):

        self.thermo = ThermoAnalyzer()

        self.conservation = PrimerConservation()

        self.snp = SNPScorer()

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
            forward["best_score"]
            + reverse["best_score"]
        ) / 2.0

        snp = (
            self.snp.score(pair.forward, snps)
            + self.snp.score(pair.reverse, snps)
        ) / 2.0

        final_score = round(
            thermo * 0.5
            + conservation * 0.3
            + snp * 0.2,
            2,
        )

        pair.score = final_score

        return {

            "thermo": round(
                thermo,
                2,
            ),

            "conservation": round(
                conservation,
                2,
            ),

            "snp": round(
                snp,
                2,
            ),

            "score": final_score,

            "final_score": final_score,

        }
