"""
Sequence manager.
"""

from pathlib import Path

from primerforge.io.cache import SequenceCache
from primerforge.io.downloader import NCBIDownloader


class SequenceManager:
    """
    Download sequences only when necessary.
    """

    def __init__(self):

        self.cache = SequenceCache()

        self.downloader = NCBIDownloader()

    def get(
        self,
        species: str,
        gene: str,
    ) -> Path:

        path = self.cache.path(
            species,
            gene,
        )

        if path.exists():
            return path

        ids = self.downloader.search(
            species,
            gene,
        )

        self.downloader.fetch(
            ids,
            path,
        )

        return path
