"""
Primer conservation tests.
"""

from pathlib import Path

from primerforge.models.primer import Primer

from primerforge.validation.conservation import PrimerConservation


def test_conservation():

    primer = Primer(
        sequence="CCCTTCATCATTGCAGCTGC",
        start=1,
        end=20,
        strand="+",
        length=20,
    )

    result = PrimerConservation().evaluate(
        primer,
        Path(
            "data/alignments/alignment.fasta"
        ),
    )

    assert result["matches"] > 0

    assert result["coverage"] > 0
