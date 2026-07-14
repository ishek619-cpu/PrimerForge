"""
Species designer tests.
"""

from envoprimer.species.designer import SpeciesPrimerDesigner


def test_designer():

    designer = SpeciesPrimerDesigner()

    assert designer.discovery is not None
