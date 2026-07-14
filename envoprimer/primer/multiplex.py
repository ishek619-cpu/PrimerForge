"""
Multiplex primer compatibility analysis.
"""

from envoprimer.models.pair import PrimerPair
from envoprimer.primer3.thermo import ThermoAnalyzer


class MultiplexAnalyzer:
    """
    Analyse compatibility among primer pairs for multiplex PCR.
    """

    def __init__(self):

        self.thermo = ThermoAnalyzer()

    def evaluate(
        self,
        pairs: list[PrimerPair],
    ):

        compatible = []
        conflicts = []

        for i in range(len(pairs)):

            for j in range(i + 1, len(pairs)):

                p1 = pairs[i]
                p2 = pairs[j]

                dimer = max(

                    self.thermo.heterodimer(
                        p1.forward,
                        p2.forward,
                    ),

                    self.thermo.heterodimer(
                        p1.forward,
                        p2.reverse,
                    ),

                    self.thermo.heterodimer(
                        p1.reverse,
                        p2.forward,
                    ),

                    self.thermo.heterodimer(
                        p1.reverse,
                        p2.reverse,
                    ),

                )

                if dimer < 10:

                    compatible.append(
                        (
                            i,
                            j,
                        )
                    )

                else:

                    conflicts.append(
                        (
                            i,
                            j,
                            round(dimer, 2),
                        )
                    )

        return {

            "compatible": compatible,

            "conflicts": conflicts,

            "compatible_pairs": len(compatible),

            "conflicting_pairs": len(conflicts),

        }
