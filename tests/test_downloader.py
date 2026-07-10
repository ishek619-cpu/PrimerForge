"""
Tests for NCBI downloader.
"""

from pathlib import Path
from unittest.mock import MagicMock

from primerforge.io.downloader import NCBIDownloader


def test_search():

    downloader = NCBIDownloader()

    downloader.search = MagicMock(
        return_value=[
            "NC_013663.1",
            "PQ810007.1",
            "PQ810734.1",
        ]
    )

    ids = downloader.search(
        "Oreochromis niloticus",
        "CYTB",
        limit=5,
    )

    assert isinstance(ids, list)

    assert len(ids) == 3

    assert ids[0] == "NC_013663.1"


def test_fetch(tmp_path: Path):

    downloader = NCBIDownloader()

    output = tmp_path / "cytb.fasta"

    def fake_fetch(ids, outfile):

        outfile.write_text(
            """>NC_013663.1 Oreochromis niloticus CYTB
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
        )

        return outfile

    downloader.fetch = MagicMock(
        side_effect=fake_fetch,
    )

    ids = [
        "NC_013663.1",
        "PQ810007.1",
    ]

    downloader.fetch(
        ids,
        output,
    )

    assert output.exists()

    text = output.read_text()

    assert text.startswith(">")

    assert "ATGCGATCGATCGATCGATC" in text
