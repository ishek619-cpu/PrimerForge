"""
GenBank annotation parser.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


@dataclass(slots=True)
class Gene:

    name: str

    feature_type: str

    start: int

    end: int

    strand: int | None

    qualifiers: dict = field(
        default_factory=dict,
    )

    @property
    def length(self) -> int:

        return self.end - self.start + 1

    def contains(
        self,
        position: int,
    ) -> bool:

        return self.start <= position <= self.end

    def to_dict(self):

        return {

            "name": self.name,

            "type": self.feature_type,

            "start": self.start,

            "end": self.end,

            "length": self.length,

            "strand": self.strand,

            "qualifiers": self.qualifiers,

        }


class AnnotationParser:
    """
    Parse GenBank annotations.
    """

    def __init__(self):

        self.genes: list[Gene] = []

    def load(
        self,
        annotation: Path,
    ) -> list[Gene]:

        self.genes.clear()

        record = SeqIO.read(
            annotation,
            "genbank",
        )

        for feature in record.features:

            if feature.type not in {

                "gene",

                "CDS",

                "rRNA",

                "tRNA",

                "misc_feature",

            }:

                continue

            self.genes.append(
                self._parse_feature(
                    feature,
                )
            )

        return self.genes

    def _parse_feature(
        self,
        feature: SeqFeature,
    ) -> Gene:

        qualifiers = dict(
            feature.qualifiers,
        )

        name = (

            qualifiers.get(
                "gene",
                qualifiers.get(
                    "product",
                    ["unknown"],
                ),
            )[0]

        )

        return Gene(

            name=name,

            feature_type=feature.type,

            start=int(
                feature.location.start
            ) + 1,

            end=int(
                feature.location.end
            ),

            strand=feature.location.strand,

            qualifiers=qualifiers,

        )

    def get_gene(
        self,
        name: str,
    ) -> Gene | None:

        for gene in self.genes:

            if gene.name.lower() == name.lower():

                return gene

        return None

    def feature_types(
        self,
    ) -> list[str]:

        return sorted(

            {

                gene.feature_type

                for gene in self.genes

            }

        )

    def gene_names(
        self,
    ) -> list[str]:

        return sorted(

            {

                gene.name

                for gene in self.genes

            }

        )
