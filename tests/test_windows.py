"""
Diagnostic window tests.
"""

from primerforge.models.region import Region
from primerforge.species.windows import DiagnosticWindow


def test_to_region():

    window = DiagnosticWindow(

        start=100,

        end=280,

        sequence="A" * 180,

        diagnostic_sites=[105, 141, 188],

    )

    region = window.to_region()

    assert isinstance(region, Region)

    assert region.start == 100

    assert region.end == 280

    assert region.length == 180

    assert region.score == 3.0
