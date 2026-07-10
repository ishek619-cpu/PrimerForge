"""
Tests for the Reference model.
"""

from pathlib import Path

import pytest

from primerforge.reference.reference import Reference


def test_reference_loading(tmp_path: Path):

    fasta = tmp_path / "reference.fasta"

    fasta.write_text(
        """>NC_000001
ATGCGTACGTAGCTAGCTAG
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
    )

    assert reference.accession == "NC_000001"

    assert reference.length == 20

    assert len(reference) == 20

    assert reference.sequence == "ATGCGTACGTAGCTAGCTAG"

    assert reference.has_alignment is False

    assert reference.has_annotation is False

    assert reference.has_metadata is False


def test_reference_alignment(tmp_path: Path):

    fasta = tmp_path / "reference.fasta"

    alignment = tmp_path / "alignment.fasta"

    fasta.write_text(
        """>ref
ATGC
""",
        encoding="utf-8",
    )

    alignment.write_text(
        """>ref
ATGC
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
        alignment=alignment,
    )

    assert reference.has_alignment

    assert not reference.has_annotation


def test_reference_summary(tmp_path: Path):

    fasta = tmp_path / "reference.fasta"

    fasta.write_text(
        """>ABC123 Example reference
ATGCATGC
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
    )

    summary = reference.summary()

    assert summary["accession"] == "ABC123"

    assert summary["length"] == 8

    assert "description" in summary

    assert "gene" in summary

    assert "species" in summary


def test_missing_reference():

    with pytest.raises(FileNotFoundError):

        Reference(
            fasta=Path("does_not_exist.fasta"),
        )


def test_repr(tmp_path: Path):

    fasta = tmp_path / "reference.fasta"

    fasta.write_text(
        """>ref
ATGC
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
    )

    text = repr(reference)

    assert "Reference" in text

    assert "length=4" in text
