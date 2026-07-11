from pathlib import Path

from primerforge.species.downloader import DownloadResult


def test_download_result():

    result = DownloadResult(

        species="Oreochromis niloticus",

        marker="COI",

        accessions=["1", "2", "3"],

        output_fasta=Path("target.fasta"),

    )

    assert result.species == "Oreochromis niloticus"

    assert result.marker == "COI"

    assert len(result.accessions) == 3

    assert result.output_fasta.name == "target.fasta"
