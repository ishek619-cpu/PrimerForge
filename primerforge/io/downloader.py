"""
NCBI sequence downloader.
"""

from pathlib import Path
from urllib.parse import quote

import requests


class NCBIDownloader:
    """
    Download nucleotide sequences from NCBI.
    """

    BASE = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/"
    )

    def search(
        self,
        species: str,
        gene: str,
        limit: int = 100,
    ) -> list[str]:

        term = quote(
            f'{species}[Organism] AND {gene}[Gene]'
        )

        url = (
            self.BASE
            + "esearch.fcgi"
            + f"?db=nucleotide&retmode=json"
            + f"&retmax={limit}"
            + f"&term={term}"
        )

        r = requests.get(
            url,
            timeout=60,
        )

        r.raise_for_status()

        return (
            r.json()["esearchresult"]["idlist"]
        )

    def fetch(
        self,
        accessions: list[str],
        output: Path,
    ) -> Path:

        ids = ",".join(accessions)

        url = (
            self.BASE
            + "efetch.fcgi"
            + "?db=nucleotide"
            + "&rettype=fasta"
            + "&retmode=text"
            + f"&id={ids}"
        )

        r = requests.get(
            url,
            timeout=120,
        )

        r.raise_for_status()

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        output.write_text(
            r.text,
        )

        return output
