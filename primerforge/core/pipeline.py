"""
PrimerForge Pipeline Engine
"""

from pathlib import Path

from primerforge.analysis.extract import GeneExtractor
from primerforge.core.config import load_config
from primerforge.core.logger import get_logger
from primerforge.io.downloader import NCBIDownloader


class Pipeline:

    def __init__(
        self,
        config_file: Path,
    ):

        self.config = load_config(config_file)

        self.logger = get_logger(__name__)

        self.downloader = NCBIDownloader(
            Path("data/genomes"),
            email="ishek619@gmail.com",
        )

        self.extractor = GeneExtractor()

    def run(self):

        self.logger.info("=" * 60)
        self.logger.info("PrimerForge Pipeline")
        self.logger.info("=" * 60)

        target = self.config.reference_taxon

        marker = self.config.marker.gene

        self.logger.info(
            f"Searching {target.name}"
        )

        ids = self.downloader.search_taxon(
            target.rank,
            target.name,
            marker,
            self.config.download.max_records,
        )

        self.logger.info(
            f"Retrieved {len(ids)} records"
        )

        accessions = self.downloader.fetch_accessions(ids)

        self.logger.info(
            f"Retrieved {len(accessions)} accessions"
        )

        for accession in accessions:

            gb = self.downloader.download_genbank(
                accession,
            )

            gene = self.extractor.extract_gene(
                gb,
                marker,
            )

            self.extractor.save_gene(
                gene,
                Path("data/genes"),
            )

        self.logger.info("Pipeline complete.")
