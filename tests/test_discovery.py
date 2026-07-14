"""
Tests for primer discovery.
"""

from pathlib import Path

from envoprimer.primer.discovery import PrimerDiscovery


def test_primer_discovery():

    pairs = PrimerDiscovery().discover(
        Path("data/genes/NC_013663_CYTB.fasta"),
    )

    assert len(pairs) > 0

    best = pairs[0]

    assert best.score > 0

    assert 80 <= best.product_size <= 250

    assert best.forward.tm > 50

    assert best.reverse.tm > 50
