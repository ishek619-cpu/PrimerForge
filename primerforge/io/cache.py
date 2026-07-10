"""
Sequence cache.
"""

from pathlib import Path


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
