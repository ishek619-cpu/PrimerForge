"""
Nearest-neighbor thermodynamic calculations.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ThermodynamicResult:
    """
    Thermodynamic properties for a duplex.
    """

    delta_h: float
    delta_s: float
    delta_g: float
    tm: float


class NearestNeighborCalculator:
    """
    Simplified nearest-neighbor thermodynamic calculator.

    This class provides the framework for SantaLucia-style
    nearest-neighbor calculations. The parameter table will
    be expanded in later sprints.
    """

    #
    # Initial perfectly matched parameters.
    # Values are placeholders and will be replaced by the
    # complete SantaLucia parameter set.
    #
    PARAMETERS = {

        "AA": (-7.9, -22.2),
        "TT": (-7.9, -22.2),

        "AT": (-7.2, -20.4),
        "TA": (-7.2, -21.3),

        "CA": (-8.5, -22.7),
        "TG": (-8.5, -22.7),

        "GT": (-8.4, -22.4),
        "AC": (-8.4, -22.4),

        "CT": (-7.8, -21.0),
        "AG": (-7.8, -21.0),

        "GA": (-8.2, -22.2),
        "TC": (-8.2, -22.2),

        "CG": (-10.6, -27.2),
        "GC": (-9.8, -24.4),

        "GG": (-8.0, -19.9),
        "CC": (-8.0, -19.9),

    }

    def calculate(
        self,
        sequence: str,
    ) -> ThermodynamicResult:

        sequence = sequence.upper()

        delta_h = 0.0
        delta_s = 0.0

        for i in range(len(sequence) - 1):

            dinucleotide = sequence[i:i + 2]

            h, s = self.PARAMETERS.get(
                dinucleotide,
                (-7.0, -20.0),
            )

            delta_h += h
            delta_s += s

        delta_g = delta_h - (310.15 * delta_s / 1000)

        tm = (
            (1000 * delta_h)
            /
            (delta_s + 1)
        )

        return ThermodynamicResult(
            delta_h=round(delta_h, 2),
            delta_s=round(delta_s, 2),
            delta_g=round(delta_g, 2),
            tm=round(tm, 2),
        )
