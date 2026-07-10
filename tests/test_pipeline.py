"""
Pipeline integration test.
"""

import pytest
from pathlib import Path

from primerforge.core.pipeline import Pipeline


@pytest.mark.slow
def test_pipeline():

    pipeline = Pipeline()

    pipeline.run(
        Path(
            "configs/example.yaml"
        )
    )

    assert Path(
        "results/alignment.fasta"
    ).exists()

    assert Path(
        "results/primers.csv"
    ).exists()

    assert Path(
        "results/primers.json"
    ).exists()

    assert Path(
        "results/index.html"
    ).exists()

    assert Path(
        "results/primers.xlsx"
    ).exists()
