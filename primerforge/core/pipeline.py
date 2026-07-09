"""
PrimerForge Pipeline Engine
"""

from pathlib import Path

from primerforge.core.config import load_config
from primerforge.core.logger import get_logger


class Pipeline:

    def __init__(
        self,
        config_file: Path,
    ):

        self.config = load_config(
            config_file,
        )

        self.logger = get_logger(__name__)

    def run(self):

        self.logger.info(
            "=" * 60
        )

        self.logger.info(
            "PrimerForge Pipeline Started"
        )

        self.logger.info(
            "=" * 60
        )

        self.logger.info(
            f"Project : {self.config.project.name}"
        )

        self.logger.info(
            f"Reference : {self.config.reference_taxon.name}"
        )

        self.logger.info(
            f"Marker : {self.config.marker.gene}"
        )

        self.logger.info(
            "Pipeline initialized successfully."
        )
