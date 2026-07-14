"""
Common species module models.

Every module in primerforge.species should use these shared
objects so the whole pipeline has one consistent API.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class SpeciesContext:
    """
    Shared state passed through the species pipeline.
    """

    species: str

    marker: str

    taxid: str = ""

    genus: str = ""

    family: str = ""

    order: str = ""

    target_accessions: list[str] = field(
        default_factory=list,
    )

    relative_species: list[str] = field(
        default_factory=list,
    )

    target_fasta: Path | None = None

    background_fasta: Path | None = None

    cleaned_target: Path | None = None

    cleaned_background: Path | None = None

    merged_fasta: Path | None = None

    alignment: Path | None = None

    diagnostics: list = field(
        default_factory=list,
    )

    metadata: dict = field(
        default_factory=dict,
    )

    @property
    def ready_for_alignment(self) -> bool:

        return (
            self.cleaned_target is not None
            and
            self.cleaned_background is not None
        )

    @property
    def ready_for_diagnostics(self) -> bool:

        return self.alignment is not None

    def summary(self) -> dict:

        return {

            "species": self.species,

            "marker": self.marker,

            "taxid": self.taxid,

            "genus": self.genus,

            "family": self.family,

            "order": self.order,

            "targets": len(
                self.target_accessions,
            ),

            "relatives": len(
                self.relative_species,
            ),

        }
