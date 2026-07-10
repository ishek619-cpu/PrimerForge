"""
Reference sequence model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO


@dataclass(slots=True)
class Reference:
    """
    Represents the reference sequence used throughout PrimerForge.

    Every major module (primer discovery, population analysis,
    specificity, PCR prediction, reporting) should receive a
    Reference object instead of separate FASTA, alignment and
    annotation paths.
    """

    ####################################################################
    # Input files
    ####################################################################

    fasta: Path

    alignment: Path | None = None

    annotation: Path | None = None

    ####################################################################
    # Metadata
    ####################################################################

    accession: str = ""

    description: str = ""

    sequence: str = ""

    length: int = 0

    species: str = ""

    gene: str = ""

    molecule_type: str = ""

    metadata: dict = field(
        default_factory=dict,
    )

    ####################################################################
    # Initialization
    ####################################################################

    def __post_init__(self):

        if not self.fasta.exists():

            raise FileNotFoundError(
                self.fasta,
            )

        record = next(
            SeqIO.parse(
                self.fasta,
                "fasta",
            )
        )

        self.accession = record.id

        self.description = record.description

        self.sequence = str(
            record.seq,
        )

        self.length = len(
            self.sequence,
        )

        if self.length == 0:

            raise ValueError(
                "Reference sequence is empty.",
            )

    ####################################################################
    # Properties
    ####################################################################

    @property
    def has_alignment(self) -> bool:

        return self.alignment is not None

    @property
    def has_annotation(self) -> bool:

        return self.annotation is not None

    @property
    def has_metadata(self) -> bool:

        return len(
            self.metadata,
        ) > 0

    ####################################################################
    # Utility methods
    ####################################################################

    def summary(self) -> dict:

        return {

            "accession": self.accession,

            "description": self.description,

            "length": self.length,

            "species": self.species,

            "gene": self.gene,

            "molecule_type": self.molecule_type,

            "alignment": self.alignment,

            "annotation": self.annotation,

        }

    def __len__(self):

        return self.length

    def __repr__(self):

        return (
            f"Reference("
            f"accession='{self.accession}', "
            f"length={self.length})"
        )
