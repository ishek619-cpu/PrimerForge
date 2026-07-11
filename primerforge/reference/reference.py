"""
Reference sequence model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

from primerforge.reference.annotation import (
    AnnotationParser,
    Gene,
)


@dataclass(slots=True)
class Reference:
    """
    Represents a biological reference sequence together with
    its alignment and annotation.
    """

    ####################################################################
    # Input files
    ####################################################################

    fasta: Path

    alignment: Path | None = None

    annotation: Path | None = None

    ####################################################################
    # Sequence information
    ####################################################################

    accession: str = ""

    description: str = ""

    sequence: str = ""

    length: int = 0

    ####################################################################
    # Biological metadata
    ####################################################################

    species: str = ""

    gene: str = ""

    molecule_type: str = ""

    metadata: dict = field(
        default_factory=dict,
    )

    ####################################################################
    # Annotation
    ####################################################################

    genes: list[Gene] = field(
        default_factory=list,
    )

    ####################################################################
    # Initialization
    ####################################################################

    def __post_init__(self):

        #
        # Validate FASTA
        #
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
                "Reference sequence is empty."
            )

        #
        # Load annotation (optional)
        #
        self.genes = []

        if self.annotation is not None:

            if self.annotation.exists():

                parser = AnnotationParser()

                self.genes = parser.load(
                    self.annotation,
                )

    ####################################################################
    # Properties
    ####################################################################

    @property
    def has_alignment(self) -> bool:

        return self.alignment is not None

    @property
    def has_annotation(self) -> bool:

        return len(
            self.genes,
        ) > 0

    @property
    def has_metadata(self) -> bool:

        return len(
            self.metadata,
        ) > 0

    ####################################################################
    # Annotation helpers
    ####################################################################

    def get_gene(
        self,
        name: str,
    ) -> Gene | None:

        for gene in self.genes:

            if gene.name.lower() == name.lower():

                return gene

        return None

    def gene_names(
        self,
    ) -> list[str]:

        return sorted(

            {

                gene.name

                for gene in self.genes

            }

        )

    ####################################################################
    # Summary
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

            "genes": len(
                self.genes,
            ),

        }

    ####################################################################
    # Convenience
    ####################################################################

    def __len__(self):

        return self.length

    def __repr__(self):

        return (

            f"Reference("

            f"accession='{self.accession}', "

            f"length={self.length}, "

            f"genes={len(self.genes)})"

        )
