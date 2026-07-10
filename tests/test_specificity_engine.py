"""
Tests for specificity engine.
"""

from pathlib import Path
from unittest.mock import MagicMock

from primerforge.models.pair import PrimerPair
from primerforge.models.primer import Primer

from primerforge.specificity.engine import SpecificityEngine
from primerforge.specificity.models import SpecificityResult


def test_specificity_engine(tmp_path: Path):

    engine = SpecificityEngine(
        database="dummy_db",
        target_species="Oreochromis niloticus",
    )

    fake = SpecificityResult(
        specificity_score=100.0,
        passed=True,
    )

    engine.analyzer.analyse = MagicMock(
        return_value=fake,
    )

    pair = PrimerPair(

        forward=Primer(
            sequence="AAAAAAAAAAAAAAAAAAAA",
            start=1,
            end=20,
            strand="+",
            length=20,
        ),

        reverse=Primer(
            sequence="TTTTTTTTTTTTTTTTTTTT",
            start=100,
            end=120,
            strand="-",
            length=20,
        ),

        product_size=120,
    )

    result = engine.evaluate_pair(
        pair,
        tmp_path,
    )

    assert result.passed

    assert pair.passed_specificity

    assert pair.specificity_score == 100.0

    assert pair.specificity_result is not None
