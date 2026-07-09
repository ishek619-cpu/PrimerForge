"""
PrimerForge Gene Extraction Engine

Extract genes from GenBank files and save them as FASTA.
"""

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from primerforge.core.logger import get_logger


class GeneExtractor:
    """
    Extract genes from GenBank files.
    """

    def __init__(self):

        self.logger = get_logger(__name__)

    def extract_gene(
        self,
        genbank_file: Path,
        gene_name: str,
    ) -> SeqRecord | None:
        """
        Extract a gene from a GenBank file.
        """

        record = SeqIO.read(
            genbank_file,
            "genbank",
        )

        for feature in record.features:

            if feature.type != "gene":
                continue

            if "gene" not in feature.qualifiers:
                continue

            if feature.qualifiers["gene"][0].upper() != gene_name.upper():
                continue

            sequence = feature.extract(record.seq)

            return SeqRecord(
                sequence,
                id=record.id,
                name=record.name,
                description=f"{record.annotations['organism']} {gene_name}",
            )

        return None

    def save_gene(
        self,
        gene_record: SeqRecord,
        output_dir: Path,
    ) -> Path:
        """
        Save extracted gene as FASTA.
        """

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        filename = (
            gene_record.description.replace(" ", "_")
            + ".fasta"
        )

        outfile = output_dir / filename

        SeqIO.write(
            gene_record,
            outfile,
            "fasta",
        )

        self.logger.info(
            f"Saved {outfile}"
        )

        return outfile
