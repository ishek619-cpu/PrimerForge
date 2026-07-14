"""
Tests for the ReferenceValidator.
"""

from pathlib import Path

from envoprimer.reference.reference import Reference
from envoprimer.reference.validator import ReferenceValidator


def test_valid_reference(tmp_path: Path):

    fasta = tmp_path / "reference.fasta"

    alignment = tmp_path / "alignment.fasta"

    fasta.write_text(
        """>ref
ATGCATGCATGC
""",
        encoding="utf-8",
    )

    alignment.write_text(
        """>ref
ATGCATGCATGC
>seq2
ATGCATGCATGC
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
        alignment=alignment,
    )

    validator = ReferenceValidator()

    result = validator.validate(
        reference,
    )

    assert result.passed

    assert result.messages == []


def test_missing_alignment(tmp_path: Path):

    fasta = tmp_path / "reference.fasta"

    fasta.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
        alignment=tmp_path / "missing.fasta",
    )

    validator = ReferenceValidator()

    result = validator.validate(
        reference,
    )

    assert not result.passed

    assert any(
        "Alignment file does not exist."
        in message
        for message in result.messages
    )


def test_alignment_longer_than_reference(tmp_path: Path):

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
ATGCATGC
>seq2
ATGCATGC
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
        alignment=alignment,
    )

    validator = ReferenceValidator()

    result = validator.validate(
        reference,
    )

    assert not result.passed

    assert any(
        "Alignment is longer than the reference sequence."
        in message
        for message in result.messages
    )


def test_missing_annotation(tmp_path: Path):

    fasta = tmp_path / "reference.fasta"

    fasta.write_text(
        """>ref
ATGCATGC
""",
        encoding="utf-8",
    )

    reference = Reference(
        fasta=fasta,
        annotation=tmp_path / "missing.gb",
    )

    validator = ReferenceValidator()

    result = validator.validate(
        reference,
    )

    assert not result.passed

    assert any(
        "Annotation file does not exist."
        in message
        for message in result.messages
    )
