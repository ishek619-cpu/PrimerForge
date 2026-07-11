"""
Reference validation framework.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio import AlignIO

from primerforge.reference.reference import Reference


@dataclass(slots=True)
class ValidationResult:
    """
    Result of validating a Reference object.
    """

    passed: bool

    messages: list[str]


class ReferenceValidator:
    """
    Validate reference, alignment and annotation consistency.
    """

    def validate(
        self,
        reference: Reference,
    ) -> ValidationResult:

        messages: list[str] = []

        #
        # FASTA exists
        #
        if not reference.fasta.exists():

            messages.append(
                "Reference FASTA does not exist."
            )

        #
        # Alignment validation
        #
        if reference.alignment is not None:

            if not reference.alignment.exists():

                messages.append(
                    "Alignment file does not exist."
                )

            else:

                aln = AlignIO.read(
                    reference.alignment,
                    "fasta",
                )

                if len(aln) == 0:

                    messages.append(
                        "Alignment contains no sequences."
                    )

                elif (
                    aln.get_alignment_length()
                    > reference.length
                ):

                    messages.append(
                        "Alignment is longer than the reference sequence."
                    )

        #
        # Annotation validation
        #
        # Annotation is optional, but if the user
        # supplies a path it must exist.
        #
        if reference.annotation is not None:

            if not reference.annotation.exists():

                messages.append(
                    "Annotation file does not exist."
                )

        return ValidationResult(

            passed=len(messages) == 0,

            messages=messages,

        )
