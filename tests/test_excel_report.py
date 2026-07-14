"""
Tests for Excel report.
"""

from openpyxl import load_workbook

from envoprimer.models.primer import Primer
from envoprimer.models.pair import PrimerPair

from envoprimer.report.excel import ExcelReport


def test_excel_report(tmp_path):

    forward = Primer(
        sequence="AAAAAAAAAAAAAAAAAAAA",
        start=1,
        end=20,
        strand="+",
        length=20,
        tm=60,
        gc=50,
    )

    reverse = Primer(
        sequence="CCCCCCCCCCCCCCCCCCCC",
        start=101,
        end=120,
        strand="-",
        length=20,
        tm=60,
        gc=60,
    )

    pair = PrimerPair(
        forward=forward,
        reverse=reverse,
        product_size=150,
        score=95,
    )

    outfile = tmp_path / "primers.xlsx"

    ExcelReport().write(
        [pair],
        outfile,
    )

    assert outfile.exists()

    workbook = load_workbook(outfile)

    sheet = workbook.active

    assert sheet["A2"].value == "AAAAAAAAAAAAAAAAAAAA"

    assert sheet["G2"].value == 150

    assert sheet["H2"].value == 95

