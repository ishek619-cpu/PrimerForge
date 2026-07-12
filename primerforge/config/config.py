"""
PrimerForge configuration loader.
"""

from pathlib import Path

import yaml


class Config:
    """
    Load PrimerForge YAML configuration.
    """

    def __init__(
        self,
        filename: Path,
    ):

        self.filename = Path(filename)

        with open(
            self.filename,
            "r",
            encoding="utf-8",
        ) as stream:

            self.data = yaml.safe_load(stream)

    @property
    def project(self):
        return self.data["project"]["name"]

    @property
    def version(self):
        return self.data["project"]["version"]

    @property
    def organism(self):
        return self.data["reference_taxon"]["name"]

    @property
    def taxonomic_rank(self):
        return self.data["reference_taxon"]["rank"]

    @property
    def gene(self):
        return self.data["marker"]["gene"]

    @property
    def primer(self):
        return self.data["primer"]

    @property
    def database(self):
        return self.data["database"]

    @property
    def download(self):
        return self.data["download"]

    @property
    def contrast_taxa(self):
        return self.data.get(
            "contrast_taxa",
            [],
        )

    #
    # NCBI credentials
    #

    @property
    def ncbi(self):
        return self.data.get(
            "ncbi",
            {},
        )

    @property
    def ncbi_email(self):
        return self.ncbi.get(
            "email",
            "",
        )

    @property
    def ncbi_api_key(self):
        return self.ncbi.get(
            "api_key",
            "",
        )

    def get(
        self,
        key,
        default=None,
    ):
        return self.data.get(
            key,
            default,
        )

    def __getitem__(
        self,
        key,
    ):
        return self.data[key]

    def __repr__(self):
        return (
            f"Config(project={self.project}, "
            f"organism={self.organism}, "
            f"gene={self.gene})"
        )
