"""
Tests for conserved region discovery.
"""

from pathlib import Path

from Bio import AlignIO

from primerforge.analysis.conservation import ConservedRegionFinder


def test_conservation(tmp_path: Path):

    fasta = tmp_path / "alignment.fasta"

    fasta.write_text(
        ">a\nAAAAAAAAAACCCCCCCCCCGGGGGGGGGG\n"
        ">b\nAAAAAAAAAACCCCCCCCCCGGGGGGGGGG\n"
        ">c\nAAAAAAAAAACCCCCCCCCCGGGGGGGGGA\n"
    )

    finder = ConservedRegionFinder()

    regions = finder.find(
        fasta,
        window=10,
        threshold=95.0,
    )

    assert len(regions) > 0

    best = finder.best(
        fasta,
        window=10,
    )

    assert best is not None

    assert best.score >= 95.0

    alignment = AlignIO.read(
        fasta,
        "fasta",
    )

    assert alignment.get_alignment_length() == 30
