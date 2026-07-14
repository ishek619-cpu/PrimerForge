"""
Tests for JSON report generation.
"""

import json

from envoprimer.models.primer import Primer
from envoprimer.models.pair import PrimerPair
from envoprimer.report.json import JSONReport


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

    #
    # Forward primer
    #
    assert (
        data[0]["forward"]["sequence"]
        == "ATGCATGCATGCATGCATGC"
    )

    assert (
        data[0]["forward"]["tm"]
        == 60.0
    )

    assert (
        data[0]["forward"]["gc"]
        == 50.0
    )

    #
    # Reverse primer
    #
    assert (
        data[0]["reverse"]["sequence"]
        == "CGATCGATCGATCGATCGAT"
    )

    assert (
        data[0]["reverse"]["tm"]
        == 60.5
    )

    assert (
        data[0]["reverse"]["gc"]
        == 55.0
    )

    #
    # General information
    #
    assert (
        data[0]["product_size"]
        == 120
    )

    assert (
        data[0]["score"]
        == 94.2
    )

    #
    # New Sprint 5 fields
    #
    assert "population" in data[0]

    assert "specificity" in data[0]

    assert "score_breakdown" in data[0]
