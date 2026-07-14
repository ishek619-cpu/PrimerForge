"""
IUPAC nucleotide matching utilities.

Supports all standard DNA ambiguity codes.
"""

from __future__ import annotations

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

    return IUPAC_CODES.get(
        base.upper(),
        IUPAC_CODES["N"],
    )


def matches(
    primer_base: str,
    sequence_base: str,
) -> bool:

    return bool(
        expand(primer_base)
        &
        expand(sequence_base)
    )


def compare(
    primer: str,
    sequence: str,
) -> tuple[int, float]:
    """
    Return

    mismatches,
    identity

    in one pass.
    """

    if len(primer) != len(sequence):
        raise ValueError(
            "Sequences must have equal length."
        )

    length = len(primer)

    if length == 0:
        return 0, 0.0

    matched = 0

    for a, b in zip(
        primer.upper(),
        sequence.upper(),
    ):

        if matches(a, b):
            matched += 1

    mismatches = length - matched

    identity = round(
        matched / length * 100,
        2,
    )

    return mismatches, identity


def mismatch_count(
    primer: str,
    sequence: str,
) -> int:

    return compare(
        primer,
        sequence,
    )[0]


def identity(
    primer: str,
    sequence: str,
) -> float:

    return compare(
        primer,
        sequence,
    )[1]
