"""
Tests for thermodynamic analysis.
"""

from primerforge.models.primer import Primer
from primerforge.primer3.thermo import ThermoAnalyzer


def test_hairpin_and_homodimer():

    primer = Primer(
        sequence="ATGCGTACGTAGCTAGCTAG",
        start=1,
        end=20,
        strand="+",
        length=20,
    )

    primer = ThermoAnalyzer().evaluate(
        primer,
    )

    assert primer.hairpin_score >= 0

    assert primer.self_dimer_score >= 0
