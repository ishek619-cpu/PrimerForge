"""
Species collector tests.
"""

from primerforge.species.collector import SpeciesDataset


def test_dataset():

    dataset = SpeciesDataset(

        species="Oreochromis niloticus",

        marker="COI",

        taxonomy=None,

        genus="Oreochromis",

        target_result=None,

        relative_species=[],

    )

    assert dataset.species == "Oreochromis niloticus"

    assert dataset.marker == "COI"

    assert dataset.genus == "Oreochromis"

    assert dataset.relative_species == []
