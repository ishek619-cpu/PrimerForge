"""
NCBI download cache.
"""

from __future__ import annotations

from pathlib import Path


class DownloadCache:
    """
    Stores downloaded FASTA files so they do not
    need to be downloaded again.
    """

    def __init__(
        self,
        cache_directory: Path,
    ):

        self.cache_directory = Path(
            cache_directory,
        )

        self.cache_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    def filename(
        self,
        species: str,
        marker: str,
    ) -> Path:

        name = (
            species.replace(" ", "_")
            + "_"
            + marker.upper()
            + ".fasta"
        )

        return self.cache_directory / name

    def exists(
        self,
        species: str,
        marker: str,
    ) -> bool:

        return self.filename(
            species,
            marker,
        ).exists()

    def path(
        self,
        species: str,
        marker: str,
    ) -> Path:

        return self.filename(
            species,
            marker,
        )
