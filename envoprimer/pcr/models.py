"""
PCR data models.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class PrimerBinding:

    sequence_id: int

    strand: str

    start: int

    end: int

    mismatches: int

    terminal_mismatches: int

    identity: float

    score: float


@dataclass(slots=True)
class Amplicon:

    sequence_id: int

    forward: PrimerBinding

    reverse: PrimerBinding

    length: int

    passed: bool
