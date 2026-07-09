"""
Tests for JSON report generation.
"""

import json

from primerforge.models.primer import Primer
from primerforge.models.pair import PrimerPair
from primerforge.report.json import JSONReport


def test_json_report(tmp_path):

    forward = Primer(
        sequence="ATGCATGCATGCATGCATGC",
        start=1,
        end=20,
        strand="+",
        length=20,
        tm=60.0,
        gc=50.0,
    )

    reverse = Primer(
        sequence="CGATCGATCGATCGATCGAT",
        start=101,
        end=120,
        strand="-",
        length=20,
        tm=60.5,
        gc=55.0,
    )

    pair = PrimerPair(
        forward=forward,
        reverse=reverse,
        product_size=120,
        score=94.2,
    )

    outfile = tmp_path / "report.json"

    JSONReport().write(
        [pair],
        outfile,
    )

    assert outfile.exists()

    data = json.loads(
        outfile.read_text()
    )

    assert len(data) == 1

    assert data[0]["forward"] == "ATGCATGCATGCATGCATGC"

    assert data[0]["product_size"] == 120

    assert data[0]["score"] == 94.2
