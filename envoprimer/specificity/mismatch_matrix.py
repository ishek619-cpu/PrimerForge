"""
Mismatch thermodynamic penalty matrix.

The values represent the relative severity of each
base mismatch during primer-template hybridization.

Higher values indicate more destabilizing mismatches.
"""

from typing import Final


#
# Perfect Watson-Crick pairs
#
MATCH_SCORE: Final = 0.0


#
# Penalty matrix
#
MISMATCH_MATRIX: Final[dict[tuple[str, str], float]] = {

    #
    # Perfect matches
    #
    ("A", "A"): MATCH_SCORE,
    ("C", "C"): MATCH_SCORE,
    ("G", "G"): MATCH_SCORE,
    ("T", "T"): MATCH_SCORE,

    #
    # Wobble pairs
    #
    ("G", "T"): 0.25,
    ("T", "G"): 0.25,

    #
    # Moderate penalties
    #
    ("A", "C"): 0.50,
    ("C", "A"): 0.50,

    ("A", "G"): 0.60,
    ("G", "A"): 0.60,

    ("C", "T"): 0.60,
    ("T", "C"): 0.60,

    #
    # Severe penalties
    #
    ("A", "T"): 0.90,
    ("T", "A"): 0.90,

    ("C", "G"): 0.90,
    ("G", "C"): 0.90,
}


def mismatch_penalty(
    primer_base: str,
    subject_base: str,
) -> float:
    """
    Return mismatch penalty.

    Unknown or ambiguous bases always receive
    the maximum penalty.
    """

    primer_base = primer_base.upper()
    subject_base = subject_base.upper()

    valid = {"A", "C", "G", "T"}

    #
    # Unknown / ambiguous bases
    #
    if (
        primer_base not in valid
        or subject_base not in valid
    ):
        return 1.0

    #
    # Perfect Watson-Crick match
    #
    if primer_base == subject_base:
        return MATCH_SCORE

    return MISMATCH_MATRIX.get(
        (primer_base, subject_base),
        1.0,
    )
