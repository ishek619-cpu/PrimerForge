"""
Tests for sequence manager.
"""

from pathlib import Path

from primerforge.io.sequence_manager import SequenceManager


class DummyDownloader:

    def __init__(self):

        self.calls = 0

    def search(
        self,
        species,
        gene,
        limit=100,
    ):

        self.calls += 1

        return ["12345"]

    def fetch(
        self,
        accessions,
        output: Path,
    ):

        output.write_text(
            ">seq\nATGC"
        )

        return output


def test_sequence_manager(tmp_path: Path):

    manager = SequenceManager()

    manager.cache.root = tmp_path

    manager.downloader = DummyDownloader()

    path1 = manager.get(
        "Oreochromis niloticus",
        "CYTB",
    )

    assert path1.exists()

    assert manager.downloader.calls == 1

    path2 = manager.get(
        "Oreochromis niloticus",
        "CYTB",
    )

    assert path2.exists()

    assert manager.downloader.calls == 1

    assert path1 == path2
