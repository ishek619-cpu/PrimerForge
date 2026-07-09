"""
Tests for HTML report generation.
"""

from pathlib import Path

from primerforge.models.pair import PrimerPair
from primerforge.models.primer import Primer
from primerforge.report.html import HTMLReport


def test_html_report(tmp_path):

    forward = Primer(
        sequence="ATGCGTACGTAGCTAGCTAG",
        start=1,
        end=20,
        strand="+",
        length=20,
        tm=60.0,
        gc=50.0,
        score=95.0,
    )

    reverse = Primer(
        sequence="CGATCGATCGTAGCTAGCAT",
        start=101,
        end=120,
        strand="-",
        length=20,
        tm=60.5,
        gc=55.0,
        score=94.0,
    )

    pair = PrimerPair(
        forward=forward,
        reverse=reverse,
        product_size=120,
        score=93.5,
    )

    outfile = tmp_path / "report.html"

    HTMLReport().write(
        [pair],
        outfile,
    )

    assert outfile.exists()

    html = outfile.read_text()

    assert "PrimerForge Report" in html

    assert "ATGCGTACGTAGCTAGCTAG" in html

    assert "CGATCGATCGTAGCTAGCAT" in html

    assert "93.50" in html
