"""
Tests for degenerate primer generation.
"""

from pathlib import Path

from envoprimer.primer.degenerate import (
    DegeneratePrimerGenerator,
)


def test_degenerate():

    primer = (
        DegeneratePrimerGenerator()
        .generate(
            Path(
                "data/alignments/alignment.fasta"
            ),
            start=100,
            length=20,
        )
    )

    assert len(primer) == 20

    assert isinstance(
        primer,
        str,
    )
