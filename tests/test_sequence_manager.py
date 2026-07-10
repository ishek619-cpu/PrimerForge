"""
Tests for SequenceManager.
"""

from pathlib import Path

from primerforge.io.sequence_manager import SequenceManager


class DummyManager(SequenceManager):

    def __init__(self):

        pass

    def get(
        self,
        organism,
        gene,
        max_records=1000,
    ):

        return Path(
            f"{organism}_{gene}.fasta"
        )


def test_sequence_manager():

    manager = DummyManager()

    result = manager.get(
        "Oreochromis niloticus",
        "CYTB",
    )

    assert result.name == "Oreochromis niloticus_CYTB.fasta"

    contrast = manager.get_contrast(
        [
            {
                "name": "Oreochromis",
            },
            {
                "name": "Sarotherodon",
            },
            {
                "name": "Coptodon",
            },
        ],
        "CYTB",
    )

    assert len(contrast) == 3

    assert contrast[0].name == "Oreochromis_CYTB.fasta"

    assert contrast[1].name == "Sarotherodon_CYTB.fasta"

    assert contrast[2].name == "Coptodon_CYTB.fasta"
