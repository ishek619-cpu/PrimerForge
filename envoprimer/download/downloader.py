"""
High-level sequence downloader.

Coordinates NCBI search, download and persistent caching.
"""

from __future__ import annotations

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from envoprimer.cache.engine import CacheEngine
from envoprimer.download.ncbi import NCBIClient
from envoprimer.download.query import QueryBuilder


class SequenceDownloader:

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
    ):

        self.ncbi = NCBIClient(
            email=email,
            api_key=api_key,
        )

        self.query = QueryBuilder()

        self.cache = CacheEngine()

    ####################################################################
    # Search
    ####################################################################

    def search(
        self,
        species: str,
        marker: str,
    ) -> list[str]:

        query = self.query.build(
            species,
            marker,
        )

        print()
        print("=" * 60)
        print("NCBI SEARCH")
        print("=" * 60)
        print(query)
        print()

        ids = self.ncbi.search(query)

        print(f"Found {len(ids)} records.")

        return ids

    ####################################################################
    # GenBank
    ####################################################################

    def download_genbank(
        self,
        species: str,
        marker: str,
    ) -> list[SeqRecord]:

        cache_dir = self.cache.downloads(
            species,
            marker,
        )

        cache_file = cache_dir / "records.gb"

        ############################################################
        # CACHE DEBUG
        ############################################################

        print()
        print("=" * 60)
        print("CACHE DEBUG")
        print("=" * 60)
        print(f"Species      : {species}")
        print(f"Marker       : {marker}")
        print(f"Cache dir    : {cache_dir}")
        print(f"Cache file   : {cache_file}")
        print(f"Exists       : {cache_file.exists()}")

        ############################################################
        # CACHE HIT
        ############################################################

        if cache_file.exists():

            print()
            print("=" * 60)
            print("USING CACHED GENBANK")
            print("=" * 60)

            records = list(
                SeqIO.parse(
                    cache_file,
                    "genbank",
                )
            )

            print(f"Loaded {len(records)} cached records.")

            return records

        ############################################################
        # CACHE MISS
        ############################################################

        print()
        print("=" * 60)
        print("CACHE MISS")
        print("=" * 60)

        ids = self.search(
            species,
            marker,
        )

        records = self.ncbi.fetch_genbank(
            ids,
        )

        print(
            f"Downloaded {len(records)} GenBank records."
        )

        cache_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        SeqIO.write(
            records,
            cache_file,
            "genbank",
        )

        print("Saved to cache.")

        return records

    ####################################################################
    # FASTA
    ####################################################################

    def download_fasta(
        self,
        species: str,
        marker: str,
    ) -> list[SeqRecord]:

        cache_dir = self.cache.downloads(
            species,
            marker,
        )

        cache_file = cache_dir / "records.fasta"

        ############################################################
        # CACHE HIT
        ############################################################

        if cache_file.exists():

            print()
            print("=" * 60)
            print("USING CACHED FASTA")
            print("=" * 60)

            records = list(
                SeqIO.parse(
                    cache_file,
                    "fasta",
                )
            )

            print(f"Loaded {len(records)} cached records.")

            return records

        ############################################################
        # CACHE MISS
        ############################################################

        ids = self.search(
            species,
            marker,
        )

        records = self.ncbi.fetch_fasta(
            ids,
        )

        print(
            f"Downloaded {len(records)} FASTA records."
        )

        cache_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        SeqIO.write(
            records,
            cache_file,
            "fasta",
        )

        print("Saved to cache.")

        return records
