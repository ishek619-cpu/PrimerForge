from pathlib import Path

from envoprimer.primer.discovery import PrimerDiscovery
from envoprimer.report.csv import CSVReport


def test_csv_report():

    pairs = PrimerDiscovery().discover(
        Path(
            "data/genes/NC_013663_CYTB.fasta"
        )
    )

    output = CSVReport().write(
        pairs,
        Path(
            "results/primers.csv"
        ),
    )

    assert output.exists()
