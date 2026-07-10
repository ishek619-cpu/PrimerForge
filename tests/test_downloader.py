"""
Tests for NCBI downloader.
"""

from pathlib import Path

from primerforge.io.downloader import (
    NCBIDownloader,
)


def test_search():

    downloader = NCBIDownloader()

    ids = downloader.search(
        "Oreochromis niloticus",
        "CYTB",
        limit=5,
    )

    assert isinstance(ids, list)

    assert len(ids) > 0


def test_fetch(tmp_path: Path):

    downloader = NCBIDownloader()

    ids = downloader.search(
        "Oreochromis niloticus",
        "CYTB",
        limit=2,
    )

    output = tmp_path / "cytb.fasta"

    downloader.fetch(
        ids,
        output,
    )

    assert output.exists()

    text = output.read_text()

    assert text.startswith(">")
