"""
PrimerForge Configuration Engine

Loads and validates YAML configuration files.
"""

from dataclasses import dataclass
from pathlib import Path

import yaml


# ---------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------

@dataclass
class ProjectConfig:
    name: str
    version: str


@dataclass
class TaxonConfig:
    rank: str
    name: str


@dataclass
class MarkerConfig:
    gene: str


@dataclass
class DatabaseConfig:
    source: str
    db: str


@dataclass
class DownloadConfig:
    include_complete_genomes: bool
    include_partial_sequences: bool
    max_records: int


@dataclass
class PrimerConfig:
    min_length: int
    max_length: int
    min_tm: float
    max_tm: float
    min_gc: float
    max_gc: float
    min_product: int
    max_product: int


@dataclass
class PrimerForgeConfig:
    project: ProjectConfig
    reference_taxon: TaxonConfig
    contrast_taxa: list[TaxonConfig]
    marker: MarkerConfig
    database: DatabaseConfig
    download: DownloadConfig
    primer: PrimerConfig


# ---------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------

def load_config(path: Path) -> PrimerForgeConfig:

    if not path.exists():
        raise FileNotFoundError(path)

    with open(path, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)

    return PrimerForgeConfig(
        project=ProjectConfig(**cfg["project"]),
        reference_taxon=TaxonConfig(**cfg["reference_taxon"]),
        contrast_taxa=[
            TaxonConfig(**taxon)
            for taxon in cfg["contrast_taxa"]
        ],
        marker=MarkerConfig(**cfg["marker"]),
        database=DatabaseConfig(**cfg["database"]),
        download=DownloadConfig(**cfg["download"]),
        primer=PrimerConfig(**cfg["primer"]),
    )
