"""
Tests for primer coverage analysis.
"""

from pathlib import Path

from envoprimer.coverage.analyzer import CoverageAnalyzer
from envoprimer.models.pair import PrimerPair
from envoprimer.models.primer import Primer


def test_coverage(tmp_path: Path):

    fasta = tmp_path / "targets.fasta"

    fasta.write_text(
        """>seq1
ATGCGTACGTAGCTAGCTAGCTAGCGTACGATCG
>seq2
TTTTATGCGTACGTAGCTAGCTAGCGGGGGGGGG
>seq3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
""",
        encoding="utf-8",
    )

    forward = Primer(
        sequence="ATGCGTAC",
        start=1,
        end=8,
        strand="+",
        length=8,
        tm=60.0,
        gc=50.0,
    )

    reverse = Primer(
        sequence="CTAGCTAG",
        start=101,
        end=108,
        strand="-",
        length=8,
        tm=60.0,
        gc=50.0,
    )

    pair = PrimerPair(
        forward=forward,
        reverse=reverse,
        product_size=150,
    )

    analyzer = CoverageAnalyzer()

    result = analyzer.analyse(
        pair,
        fasta,
    )

    assert result["total"] == 3

    assert result["forward_count"] == 2
    assert result["reverse_count"] == 2
    assert result["pair_count"] == 2

    assert result["forward"] == 66.67
    assert result["reverse"] == 66.67
    assert result["pair"] == 66.67
