"""
Data models for the EnvoPrimer v2 PCR engine.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class BindingSite:
    """
    Represents one primer binding site on a sequence.
    """

    sequence_id: int
    strand: str              # "+" or "-"
    start: int               # 0-based inclusive
    end: int                 # 0-based exclusive

    identity: float
    mismatches: int
    terminal_mismatches: int

    score: float

    primer_alignment: str
    template_alignment: str


@dataclass(slots=True)
class PCRProduct:
    """
    Represents one predicted PCR amplicon.
    """

    sequence_id: int

    forward_site: BindingSite
    reverse_site: BindingSite

    start: int
    end: int
    size: int

    predicted: bool
