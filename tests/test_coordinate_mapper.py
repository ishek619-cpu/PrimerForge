from pathlib import Path

from primerforge.alignment.coordinate_mapper import (
    CoordinateMapper,
)


def test_reference_to_alignment(tmp_path: Path):

    fasta = tmp_path / "alignment.fasta"

    fasta.write_text(
        """>reference
ATGC--ATGC
>sample1
ATGCTTATGC
""",
        encoding="utf-8",
    )

    mapper = CoordinateMapper(fasta)

    assert mapper.reference_to_alignment(0) == 0
    assert mapper.reference_to_alignment(1) == 1
    assert mapper.reference_to_alignment(2) == 2
    assert mapper.reference_to_alignment(3) == 3

    #
    # Two gaps appear after reference position 3.
    #
    assert mapper.reference_to_alignment(4) == 6
    assert mapper.reference_to_alignment(5) == 7
    assert mapper.reference_to_alignment(6) == 8
    assert mapper.reference_to_alignment(7) == 9


def test_alignment_to_reference(tmp_path: Path):

    fasta = tmp_path / "alignment.fasta"

    fasta.write_text(
        """>reference
ATGC--ATGC
>sample1
ATGCTTATGC
""",
        encoding="utf-8",
    )

    mapper = CoordinateMapper(fasta)

    assert mapper.alignment_to_reference(0) == 0
    assert mapper.alignment_to_reference(3) == 3

    #
    # Gap columns still map to the previous
    # ungapped reference coordinate.
    #
    assert mapper.alignment_to_reference(4) == 4
    assert mapper.alignment_to_reference(5) == 4

    assert mapper.alignment_to_reference(6) == 4
    assert mapper.alignment_to_reference(7) == 5
    assert mapper.alignment_to_reference(8) == 6
    assert mapper.alignment_to_reference(9) == 7
