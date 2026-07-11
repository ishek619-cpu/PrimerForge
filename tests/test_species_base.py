"""
Tests for shared species models.
"""

from primerforge.species.base import SpeciesContext


def test_species_context():

    context = SpeciesContext(

        species="Oreochromis niloticus",

        marker="COI",

    )

    assert context.species == "Oreochromis niloticus"

    assert context.marker == "COI"

    assert context.taxid == ""

    assert context.relative_species == []

    assert context.target_accessions == []

    assert not context.ready_for_alignment

    assert not context.ready_for_diagnostics

    summary = context.summary()

    assert summary["species"] == "Oreochromis niloticus"

    assert summary["marker"] == "COI"
