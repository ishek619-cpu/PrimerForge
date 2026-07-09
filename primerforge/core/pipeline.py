"""
PrimerForge Pipeline Engine
"""

from pathlib import Path

from primerforge.analysis.align import AlignmentEngine
from primerforge.analysis.conservation import ConservationAnalyzer
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

        self.aligner = AlignmentEngine()

        self.conservation = ConservationAnalyzer()

    def run(self):

        self.logger.info("=" * 60)
        self.logger.info("PrimerForge Pipeline")
        self.logger.info("=" * 60)

        target = self.config.reference_taxon

        marker = self.config.marker.gene

        self.logger.info(f"Searching {target.name}")

        ids = self.downloader.search_taxon(
            target.rank,
            target.name,
            marker,
            self.config.download.max_records,
        )

        self.logger.info(f"Retrieved {len(ids)} records")

        accessions = self.downloader.fetch_accessions(ids)

        self.logger.info(
            f"Retrieved {len(accessions)} accessions"
        )

        genes_dir = Path("data/genes")
        genomes_dir = Path("data/genomes")
        alignments_dir = Path("data/alignments")

        genes_dir.mkdir(parents=True, exist_ok=True)
        genomes_dir.mkdir(parents=True, exist_ok=True)
        alignments_dir.mkdir(parents=True, exist_ok=True)

        for accession in accessions:

            gb = self.downloader.download_genbank(
                accession,
            )

            gene = self.extractor.extract_gene(
                gb,
                marker,
            )

            if gene is None:
                self.logger.warning(
                    f"{accession}: {marker} not found"
                )
                continue

            self.extractor.save_gene(
                gene,
                genes_dir,
            )

        self.logger.info("Combining FASTA files")

        combined = self.aligner.combine_fastas(
            genes_dir,
            alignments_dir / "all_sequences.fasta",
        )

        self.logger.info("Running MUSCLE")

        alignment = self.aligner.run_muscle(
            combined,
            alignments_dir / "alignment.fasta",
        )

        self.logger.info("Calculating conservation")

        self.conservation.summary(
            alignment,
        )

        self.logger.info("=" * 60)
        self.logger.info("Pipeline completed successfully")
        self.logger.info("=" * 60)
