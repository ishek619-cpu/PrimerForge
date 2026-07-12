"""
Alignment builder tests.
"""

from pathlib import Path

from primerforge.species.alignment import AlignmentResult


def test_alignment_result():

    result = AlignmentResult(

        target_alignment=Path("target_alignment.fasta"),

        background_alignment=Path("background_alignment.fasta"),

    )

    assert result.target_alignment.name == "target_alignment.fasta"

    assert result.background_alignment.name == "background_alignment.fasta"
