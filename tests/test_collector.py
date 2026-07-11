"""
Species collector tests.
"""

from pathlib import Path

from primerforge.species.collector import SpeciesDataset


def test_dataset():

    dataset = SpeciesDataset(

        species="Oreochromis niloticus",

        marker="COI",

        target_fasta=Path("target.fasta"),

        relative_fastas=[],

        cleaned_target=Path("cleaned.fasta"),

        cleaned_relatives=[],

        alignment=Path("alignment.fasta"),

        diagnostics=[],

    )

    assert dataset.species == "Oreochromis niloticus"

    assert dataset.marker == "COI"
