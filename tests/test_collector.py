"""
Species collector tests.
"""

from pathlib import Path

from envoprimer.species.collector import SpeciesDataset


def test_dataset():

    dataset = SpeciesDataset(

        species="Oreochromis niloticus",

        marker="COI",

        taxonomy=None,

        genus="Oreochromis",

        target_result=None,

        background_results=[],

        target_fasta=Path("target.fasta"),

        background_fasta=Path("background.fasta"),

    )

    assert dataset.species == "Oreochromis niloticus"

    assert dataset.marker == "COI"

    assert dataset.genus == "Oreochromis"

    assert dataset.background_results == []

    assert dataset.target_fasta.name == "target.fasta"

    assert dataset.background_fasta.name == "background.fasta"
