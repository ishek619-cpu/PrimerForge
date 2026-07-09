"""
Primer thermodynamic analysis using Primer3.
"""

import primer3

from primerforge.models.primer import Primer


class ThermoAnalyzer:
    """
    Calculate primer secondary structure.
    """

    @staticmethod
    def _dg(result):

        value = getattr(result, "dg", 0.0)

        if value is None:
            return 0.0

        return round(abs(value) / 1000.0, 2)

    def evaluate(
        self,
        primer: Primer,
    ) -> Primer:

        hairpin = primer3.calc_hairpin(
            primer.sequence,
        )

        homodimer = primer3.calc_homodimer(
            primer.sequence,
        )

        primer.hairpin_score = self._dg(
            hairpin,
        )

        primer.self_dimer_score = self._dg(
            homodimer,
        )

        return primer

    def heterodimer(
        self,
        forward: Primer,
        reverse: Primer,
    ) -> float:

        result = primer3.calc_heterodimer(
            forward.sequence,
            reverse.sequence,
        )

        return self._dg(result)
