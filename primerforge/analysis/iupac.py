"""
IUPAC nucleotide matching utilities.

Supports all standard DNA ambiguity codes.
"""

from __future__ import annotations

# Complete IUPAC nucleotide code table
IUPAC_CODES: dict[str, frozenset[str]] = {
    "A": frozenset({"A"}),
    "C": frozenset({"C"}),
    "G": frozenset({"G"}),
    "T": frozenset({"T"}),

    "R": frozenset({"A", "G"}),
    "Y": frozenset({"C", "T"}),
    "S": frozenset({"G", "C"}),
    "W": frozenset({"A", "T"}),
    "K": frozenset({"G", "T"}),
    "M": frozenset({"A", "C"}),

    "B": frozenset({"C", "G", "T"}),
    "D": frozenset({"A", "G", "T"}),
    "H": frozenset({"A", "C", "T"}),
    "V": frozenset({"A", "C", "G"}),

    "N": frozenset({"A", "C", "G", "T"}),

    "-": frozenset(),
}


def expand(base: str) -> frozenset[str]:
    """
    Return all nucleotides represented by an IUPAC code.

    Unknown symbols are treated as N.
    """

    return IUPAC_CODES.get(
        base.upper(),
        IUPAC_CODES["N"],
    )


def matches(
    primer_base: str,
    sequence_base: str,
) -> bool:
    """
    True if two IUPAC symbols are compatible.
    """

    return bool(
        expand(primer_base)
        &
        expand(sequence_base)
    )


def mismatch(
    primer_base: str,
    sequence_base: str,
) -> bool:
    """
    True if two bases are incompatible.
    """

    return not matches(
        primer_base,
        sequence_base,
    )


def identity(
    primer: str,
    sequence: str,
) -> float:
    """
    Percentage identity between two sequences using
    IUPAC-aware comparison.
    """

    if len(primer) != len(sequence):
        raise ValueError(
            "Sequences must have equal length."
        )

    if not primer:
        return 0.0

    matches_count = sum(
        matches(a, b)
        for a, b in zip(
            primer.upper(),
            sequence.upper(),
        )
    )

    return round(
        matches_count
        / len(primer)
        * 100.0,
        2,
    )


def mismatch_count(
    primer: str,
    sequence: str,
) -> int:
    """
    Count incompatible positions.
    """

    if len(primer) != len(sequence):
        raise ValueError(
            "Sequences must have equal length."
        )

    return sum(
        mismatch(a, b)
        for a, b in zip(
            primer.upper(),
            sequence.upper(),
        )
    )
