"""
Tests for primer specificity analysis.
"""

from pathlib import Path
from unittest.mock import MagicMock

from primerforge.models.primer import Primer
from primerforge.specificity.analyzer import SpecificityAnalyzer
from primerforge.specificity.models import (
    BlastHit,
    SpecificityResult,
)


def test_specificity_analyzer(tmp_path: Path):

    analyzer = SpecificityAnalyzer()

    analyzer.blast.search = MagicMock()

    analyzer.parser.parse = MagicMock(
        return_value=[
            BlastHit(
                accession="NC_013663.1",
                species="Oreochromis niloticus",
                identity=100.0,
                coverage=100.0,
                alignment_length=20,
                mismatches=0,
                gap_opens=0,
                qstart=1,
                qend=20,
                sstart=100,
                send=119,
                strand="plus",
                bitscore=40.0,
                evalue=1e-20,
            )
        ]
    )

    primer = Primer(
        sequence="ATGCGATCGATCGATCGATC",
        start=1,
        end=20,
        strand="+",
        length=20,
    )

    result = analyzer.analyse(
        primer=primer,
        database="dummy_db",
        blast_output=tmp_path / "blast.tsv",
        target_species="Oreochromis niloticus",
    )

    assert isinstance(result, SpecificityResult)
    assert result.passed
    assert result.specificity_score == 100.0
    assert len(result.target_hits) == 1
    assert len(result.off_target_hits) == 0
