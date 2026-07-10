"""
Sequence management for PrimerForge.
"""

from pathlib import Path

from primerforge.io.cache import SequenceCache
from primerforge.io.downloader import NCBIDownloader


class SequenceManager:
    """
    Manage downloading and caching of target sequences.
    """

    def __init__(self):

        self.cache = SequenceCache()

        self.downloader = NCBIDownloader()

    def get(
        self,
        organism: str,
        gene: str,
        max_records: int = 1000,
    ) -> Path:

        cached = self.cache.find(
            organism,
            gene,
        )

        if cached is not None:

            return cached

        print(
            f"Downloading {organism} {gene}..."
        )

        accessions = self.downloader.search(
            species=organism,
            gene=gene,
            limit=max_records,
        )

        output = self.cache.path(
            organism,
            gene,
        )

        fasta = self.downloader.fetch(
            accessions,
            output,
        )

        self.cache.store(
            organism,
            gene,
            fasta,
        )

        return self.cache.find(
            organism,
            gene,
        )

    def get_contrast(
        self,
        taxa: list,
        gene: str,
        max_records: int = 1000,
    ) -> list[Path]:

        files = []

        for taxon in taxa:

            organism = taxon["name"]

            files.append(
                self.get(
                    organism,
                    gene,
                    max_records=max_records,
                )
            )

        return files
