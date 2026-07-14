"""
Tests for primer conservation.
"""

from pathlib import Path

from envoprimer.models.primer import Primer
from envoprimer.validation.conservation import (
    PrimerConservation,
)


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
            "data/alignments/alignment.fasta",
        ),
    )

    assert "best_start" in result

    assert "window_length" in result

    assert "best_score" in result

    assert "profile" in result

    assert result["window_length"] == 20

    assert isinstance(
        result["profile"],
        list,
    )

    assert len(result["profile"]) > 0

    assert result["best_score"] >= 0.0
