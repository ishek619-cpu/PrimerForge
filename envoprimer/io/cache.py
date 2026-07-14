"""
Sequence cache.
"""

from pathlib import Path
import shutil


class SequenceCache:
    """
    Manage cached sequence files.
    """

    def __init__(
        self,
        root: Path = Path.home() / ".primerforge",
    ):

        self.root = root

        self.root.mkdir(
            parents=True,
            exist_ok=True,
        )

    def filename(
        self,
        species: str,
        gene: str,
    ) -> Path:

        safe = (
            species.lower()
            .replace(" ", "_")
            .replace("/", "_")
        )

        return self.root / f"{safe}_{gene.lower()}.fasta"

    def exists(
        self,
        species: str,
        gene: str,
    ) -> bool:

        return self.filename(
            species,
            gene,
        ).exists()

    def path(
        self,
        species: str,
        gene: str,
    ) -> Path:

        return self.filename(
            species,
            gene,
        )

    def find(
        self,
        species: str,
        gene: str,
    ) -> Path | None:

        cache_file = self.path(
            species,
            gene,
        )

        if cache_file.exists():
            return cache_file

        return None

    def store(
        self,
        species: str,
        gene: str,
        fasta: Path,
    ) -> Path:
        """
        Store a FASTA file in the cache.

        If the file is already at the cache location,
        nothing is copied.
        """

        destination = self.path(
            species,
            gene,
        )

        fasta = Path(fasta)

        if fasta.resolve() != destination.resolve():

            shutil.copy2(
                fasta,
                destination,
            )

        return destination
