"""
Tests for primer conservation.
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

    assert result["coverage"] >= 0

    assert result["exact"] >= 0

    assert result["one_mismatch"] >= 0

    assert result["two_mismatch"] >= 0

    assert result["mean_mismatches"] >= 0
