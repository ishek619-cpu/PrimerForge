"""
Reference-related exceptions.
"""

from __future__ import annotations


class ReferenceError(Exception):
    """
    Base exception for all Reference-related errors.
    """


class ReferenceFileNotFoundError(ReferenceError):
    """
    Raised when a required reference file cannot be found.
    """

    def __init__(self, path):

        super().__init__(
            f"Reference file not found: {path}"
        )


class AlignmentFileNotFoundError(ReferenceError):
    """
    Raised when an alignment file cannot be found.
    """

    def __init__(self, path):

        super().__init__(
            f"Alignment file not found: {path}"
        )


class AnnotationFileNotFoundError(ReferenceError):
    """
    Raised when an annotation file cannot be found.
    """

    def __init__(self, path):

        super().__init__(
            f"Annotation file not found: {path}"
        )


class EmptyReferenceError(ReferenceError):
    """
    Raised when the reference sequence is empty.
    """

    def __init__(self):

        super().__init__(
            "Reference sequence is empty."
        )


class CoordinateSystemError(ReferenceError):
    """
    Raised when reference and alignment use
    incompatible coordinate systems.
    """

    def __init__(
        self,
        reference_length: int,
        alignment_length: int,
    ):

        super().__init__(
            "Reference and alignment appear to use "
            "different coordinate systems "
            f"(reference={reference_length}, "
            f"alignment={alignment_length})."
        )
