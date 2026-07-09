"""
PrimerForge NCBI Downloader

Downloads sequence data and metadata from NCBI.
"""

from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez

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

            accession = record["Caption"]

            accessions.append(accession)

        self.logger.info(f"Retrieved {len(accessions)} accession numbers")

        return accessions

    def download_genbank(
        self,
        accession: str,
    ):
        raise NotImplementedError

    def download_fasta(
        self,
        accession: str,
    ):
        raise NotImplementedError

    def save_metadata(
        self,
        metadata: SequenceMetadata,
    ):
        raise NotImplementedError
