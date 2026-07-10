"""
Tests for Reference exceptions.
"""

from primerforge.reference.exceptions import (
    AlignmentFileNotFoundError,
    AnnotationFileNotFoundError,
    CoordinateSystemError,
    EmptyReferenceError,
    ReferenceError,
    ReferenceFileNotFoundError,
)


def test_reference_file_not_found():

    error = ReferenceFileNotFoundError(
        "reference.fasta",
    )

    assert isinstance(
        error,
        ReferenceError,
    )

    assert (
        str(error)
        == "Reference file not found: reference.fasta"
    )


def test_alignment_file_not_found():

    error = AlignmentFileNotFoundError(
        "alignment.fasta",
    )

    assert isinstance(
        error,
        ReferenceError,
    )

    assert (
        str(error)
        == "Alignment file not found: alignment.fasta"
    )


def test_annotation_file_not_found():

    error = AnnotationFileNotFoundError(
        "annotation.gbk",
    )

    assert isinstance(
        error,
        ReferenceError,
    )

    assert (
        str(error)
        == "Annotation file not found: annotation.gbk"
    )


def test_empty_reference():

    error = EmptyReferenceError()

    assert isinstance(
        error,
        ReferenceError,
    )

    assert (
        str(error)
        == "Reference sequence is empty."
    )


def test_coordinate_system_error():

    error = CoordinateSystemError(
        reference_length=16569,
        alignment_length=1140,
    )

    assert isinstance(
        error,
        ReferenceError,
    )

    assert (
        "different coordinate systems"
        in str(error)
    )

    assert (
        "16569"
        in str(error)
    )

    assert (
        "1140"
        in str(error)
    )
