"""
Species collection workflow.

Creates the complete species-specific dataset required for
PrimerForge.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from primerforge.species.downloader import (
    SpeciesDownloader,
)

from primerforge.species.taxonomy import (
    TaxonomyResolver,
)

from primerforge.species.relatives import (
    RelativeSpeciesFinder,
)


@dataclass(slots=True)
class SpeciesDataset:

    species: str

    marker: str

    taxonomy: object

    genus: str

    target_result: object

    relative_species: list[str]


class SpeciesCollector:

    def __init__(

        self,

        email: str,

        api_key: str | None = None,

    ):

        self.taxonomy = TaxonomyResolver(

            email=email,

            api_key=api_key,

        )

        self.relatives = RelativeSpeciesFinder(

            email=email,

            api_key=api_key,

        )

        self.downloader = SpeciesDownloader(

            email=email,

            api_key=api_key,

        )

    def collect(

        self,

        species: str,

        marker: str,

        output_directory: Path,

    ) -> SpeciesDataset:

        taxonomy = self.taxonomy.resolve(

            species,

        )

        genus = self.relatives.genus(

            taxonomy,

        )

        relative_ids = self.relatives.search(

            genus,

        )

        target = self.downloader.download(

            species=species,

            marker=marker,

            output_directory=output_directory,

        )

        return SpeciesDataset(

            species=species,

            marker=marker,

            taxonomy=taxonomy,

            genus=genus,

            target_result=target,

            relative_species=relative_ids,

        )
