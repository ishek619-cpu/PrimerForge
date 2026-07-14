"""
Tests for sequence cache.
"""

from pathlib import Path

from envoprimer.io.cache import SequenceCache


def test_cache(tmp_path: Path):

    cache = SequenceCache(tmp_path)

    path = cache.path(
        "Oreochromis niloticus",
        "CYTB",
    )

    assert path.name == "oreochromis_niloticus_cytb.fasta"

    assert cache.exists(
        "Oreochromis niloticus",
        "CYTB",
    ) is False

    path.write_text(">seq\nATGC")

    assert cache.exists(
        "Oreochromis niloticus",
        "CYTB",
    ) is True
