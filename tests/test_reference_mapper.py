"""
Tests for the universal coordinate mapper.
"""

from pathlib import Path

import pytest

from envoprimer.reference.coordinates import Coordinate
from envoprimer.reference.exceptions import CoordinateSystemError
from envoprimer.reference.mapper import CoordinateMapper
from envoprimer.reference.reference import Reference


def test_reference_to_alignment(tmp_path: Path):

    reference = tmp_path / "reference.fasta"

    alignment = tmp_path / "alignment.fasta"

    reference.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    alignment.write_text(
        """>ref
ATGCATGC
>seq2
ATGCATGC
""",
        encoding="utf-8",
    )

    ref = Reference(
        fasta=reference,
        alignment=alignment,
    )

    mapper = CoordinateMapper(
        ref,
    )

    coordinate = Coordinate(
        start=2,
        end=5,
        system="reference",
    )

    mapped = mapper.reference_to_alignment(
        coordinate,
    )

    assert mapped.start == 2

    assert mapped.end == 5

    assert mapped.system == "alignment"


def test_alignment_to_reference(tmp_path: Path):

    reference = tmp_path / "reference.fasta"

    alignment = tmp_path / "alignment.fasta"

    reference.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    alignment.write_text(
        """>ref
ATGCATGC
>seq2
ATGCATGC
""",
        encoding="utf-8",
    )

    ref = Reference(
        fasta=reference,
        alignment=alignment,
    )

    mapper = CoordinateMapper(
        ref,
    )

    coordinate = Coordinate(
        start=3,
        end=7,
        system="alignment",
    )

    mapped = mapper.alignment_to_reference(
        coordinate,
    )

    assert mapped.start == 3

    assert mapped.end == 7

    assert mapped.system == "reference"


def test_wrong_reference_system(tmp_path: Path):

    reference = tmp_path / "reference.fasta"

    alignment = tmp_path / "alignment.fasta"

    reference.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    alignment.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    ref = Reference(
        fasta=reference,
        alignment=alignment,
    )

    mapper = CoordinateMapper(
        ref,
    )

    with pytest.raises(
        CoordinateSystemError,
    ):

        mapper.reference_to_alignment(

            Coordinate(
                start=1,
                end=4,
                system="alignment",
            )

        )


def test_wrong_alignment_system(tmp_path: Path):

    reference = tmp_path / "reference.fasta"

    alignment = tmp_path / "alignment.fasta"

    reference.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    alignment.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    ref = Reference(
        fasta=reference,
        alignment=alignment,
    )

    mapper = CoordinateMapper(
        ref,
    )

    with pytest.raises(
        CoordinateSystemError,
    ):

        mapper.alignment_to_reference(

            Coordinate(
                start=1,
                end=4,
                system="reference",
            )

        )


def test_reference_without_alignment(tmp_path: Path):

    reference = tmp_path / "reference.fasta"

    reference.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    ref = Reference(
        fasta=reference,
    )

    with pytest.raises(
        ValueError,
    ):

        CoordinateMapper(
            ref,
        )
