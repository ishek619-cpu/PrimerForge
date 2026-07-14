"""
Integration test for the species workflow.
"""

from envoprimer.species.designer import SpeciesPrimerDesigner


def test_species_primer_designer_creation():

    designer = SpeciesPrimerDesigner()

    assert designer.discovery is not None
