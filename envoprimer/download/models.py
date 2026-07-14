"""
Download models used throughout EnvoPrimer.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from envoprimer.species.metadata import SequenceMetadata


@dataclass(slots=True)
class DownloadResult:
    """
    Result produced by the EnvoPrimer download pipeline.
    """

    species: str

    marker: str

    accessions: list[str]

    metadata: list[SequenceMetadata]

    output_fasta: Path
