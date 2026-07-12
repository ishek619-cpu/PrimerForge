"""
Species designer tests.
"""

from primerforge.species.designer import SpeciesPrimerDesigner


def test_designer():

    designer = SpeciesPrimerDesigner()

    assert designer.discovery is not None
