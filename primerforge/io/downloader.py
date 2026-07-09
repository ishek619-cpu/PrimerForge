"""
PrimerForge NCBI Downloader
"""

from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez
from Bio import SeqIO

from primerforge.core.logger import get_logger


@dataclass
class SequenceMetadata:
    accession: str
    organism: str
    gene: str
    length: int
    topology: str
    molecule: str


class NCBIDownloader:

    def __init__(
        self,
        output_dir: Path,
        email: str,
    ):

        self.output_dir = output_dir
        self.output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        Entrez.email = email

        self.logger = get_logger(__name__)

    def search_taxon(
        self,
        rank: str,
        name: str,
        gene: str,
        max_records: int = 100,
    ) -> list[str]:

        query = f'"{name}"[Organism] AND {gene}[Gene]'

        self.logger.info(f"Searching NCBI: {query}")

        handle = Entrez.esearch(
            db="nucleotide",
            term=query,
            retmax=max_records,
        )

        results = Entrez.read(handle)

        handle.close()

        ids = results["IdList"]

        self.logger.info(f"Found {len(ids)} records")

        return ids

    def fetch_accessions(
        self,
        ids: list[str],
    ) -> list[str]:

        self.logger.info("Fetching accession numbers")

        handle = Entrez.esummary(
            db="nucleotide",
            id=",".join(ids),
        )

        records = Entrez.read(handle)

        handle.close()

        accessions = []

        for record in records:
            accessions.append(record["Caption"])

        self.logger.info(
            f"Retrieved {len(accessions)} accession numbers"
        )

        return accessions

    def download_genbank(
        self,
        accession: str,
    ) -> Path:

        self.logger.info(f"Downloading GenBank: {accession}")

        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text",
        )

        text = handle.read()

        handle.close()

        outfile = self.output_dir / f"{accession}.gb"

        with open(outfile, "w", encoding="utf-8") as out:
            out.write(text)

        self.logger.info(f"Saved {outfile}")

        return outfile

    def download_fasta(
        self,
        accession: str,
    ) -> Path:

        self.logger.info(f"Downloading FASTA: {accession}")

        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta",
            retmode="text",
        )

        text = handle.read()

        handle.close()

        outfile = self.output_dir / f"{accession}.fasta"

        with open(outfile, "w", encoding="utf-8") as out:
            out.write(text)

        self.logger.info(f"Saved {outfile}")

        return outfile

    def extract_metadata(
        self,
        genbank_file: Path,
    ) -> SequenceMetadata:

        record = SeqIO.read(genbank_file, "genbank")

        accession = record.id

        organism = record.annotations.get(
            "organism",
            "Unknown",
        )

        length = len(record.seq)

        topology = record.annotations.get(
            "topology",
            "Unknown",
        )

        molecule = record.annotations.get(
            "molecule_type",
            "Unknown",
        )

        gene = "Unknown"

        for feature in record.features:

            if feature.type == "gene":

                if "gene" in feature.qualifiers:

                    gene = feature.qualifiers["gene"][0]

                    break

        return SequenceMetadata(
            accession=accession,
            organism=organism,
            gene=gene,
            length=length,
            topology=topology,
            molecule=molecule,
        )
