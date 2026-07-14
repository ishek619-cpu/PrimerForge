"""
Gene extraction from GenBank records.

EnvoPrimer extracts the requested marker directly from
GenBank FEATURES rather than relying on FASTA titles.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


MARKER_SYNONYMS = {

    "COI": {
        "COI",
        "CO1",
        "COXI",
        "COX1",
        "CYTOCHROME OXIDASE I",
        "CYTOCHROME C OXIDASE SUBUNIT I",
    },

    "CYTB": {
        "CYTB",
        "CYTOCHROME B",
    },

    "12S": {
        "12S",
        "12S RRNA",
        "SMALL SUBUNIT RRNA",
    },

    "16S": {
        "16S",
        "16S RRNA",
        "LARGE SUBUNIT RRNA",
    },

    "18S": {
        "18S",
        "18S RRNA",
    },

    "ITS": {
        "ITS",
        "ITS1",
        "ITS2",
        "INTERNAL TRANSCRIBED SPACER",
    },

}


@dataclass(slots=True)
class GeneRecord:

    accession: str

    gene: str

    start: int

    end: int

    strand: int

    sequence: str

    description: str


class GeneExtractor:

    def __init__(self):

        pass

    def _normalize(
        self,
        text: str,
    ) -> str:

        return (
            text.upper()
            .replace("-", " ")
            .replace("_", " ")
            .strip()
        )

    def _matches(
        self,
        requested: str,
        candidate: str,
    ) -> bool:

        requested = requested.upper()

        candidate = self._normalize(
            candidate,
        )

        if requested not in MARKER_SYNONYMS:
            return False

        return (
            candidate
            in MARKER_SYNONYMS[requested]
        )

    def read(
        self,
        genbank: Path,
    ) -> list[SeqRecord]:

        return list(
            SeqIO.parse(
                genbank,
                "genbank",
            )
        )

    def extract(
        self,
        genbank: Path,
        marker: str,
    ) -> list[GeneRecord]:

        marker = marker.upper()

        output: list[GeneRecord] = []

        for record in self.read(genbank):

            output.extend(
                self.extract_record(
                    record,
                    marker,
                )
            )

        return output

    def extract_record(
        self,
        record: SeqRecord,
        marker: str,
    ) -> list[GeneRecord]:

        genes: list[GeneRecord] = []

        for feature in record.features:

            if feature.type not in {
                "CDS",
                "gene",
                "rRNA",
            }:
                continue

            names = []

            for qualifier in (
                "gene",
                "product",
                "label",
                "standard_name",
            ):

                names.extend(
                    feature.qualifiers.get(
                        qualifier,
                        [],
                    )
                )

            if not any(
                self._matches(
                    marker,
                    name,
                )
                for name in names
            ):
                continue

            sequence = str(
                feature.extract(
                    record.seq,
                )
            )

            genes.append(
                GeneRecord(
                    accession=record.id,
                    gene=marker,
                    start=int(feature.location.start),
                    end=int(feature.location.end),
                    strand=feature.location.strand,
                    sequence=sequence,
                    description=record.description,
                )
            )

        return genes

    def write_fasta(
        self,
        genes: list[GeneRecord],
        output: Path,
    ) -> None:

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with output.open(
            "w",
            encoding="utf-8",
        ) as handle:

            for gene in genes:

                handle.write(
                    f">{gene.accession} {gene.gene}\n"
                )

                sequence = gene.sequence

                for i in range(
                    0,
                    len(sequence),
                    80,
                ):

                    handle.write(
                        sequence[i:i + 80] + "\n"
                    )
