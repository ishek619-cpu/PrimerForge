"""
Tests for the primer validation pipeline.
"""

from pathlib import Path

from envoprimer.models.pair import PrimerPair
from envoprimer.models.primer import Primer
from envoprimer.pipeline.validator import PrimerValidationPipeline


def test_pipeline_validator(tmp_path: Path):

    fasta = tmp_path / "alignment.fasta"

    fasta.write_text(
        """>seq1
ATGCGTACGTAGCTAGCTAGCTAGCGTACGATCG
>seq2
ATGCGTACGTAGCTAGCTAGCTAGCGTACGATCG
>seq3
ATGCGTACGTAGCTAGCTAGCTAGCGTACGATCG
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

    pipeline = PrimerValidationPipeline()

    result = pipeline.validate(
        pairs=[pair],
        reference_fasta=fasta,
        alignment_fasta=fasta,
        specificity_engine=None,
    )

    assert len(result) == 1

    validated = result[0]

    assert hasattr(validated, "coverage")

    assert validated.coverage["pair"] == 100.0

    assert validated.score >= 0
