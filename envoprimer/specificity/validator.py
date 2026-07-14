"""
Primer specificity validation.
"""

from pathlib import Path
from tempfile import TemporaryDirectory

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from envoprimer.specificity.blast import BlastRunner
from envoprimer.specificity.parser import BlastParser
from envoprimer.specificity.scorer import SpecificityScorer


class SpecificityValidator:
    """
    Validate primer specificity using BLAST.
    """

    def __init__(self):

        self.runner = BlastRunner()

        self.parser = BlastParser()

        self.scorer = SpecificityScorer()

    def blast_hits(
        self,
        sequence: str,
        database: str,
    ):

        with TemporaryDirectory() as tmp:

            tmp = Path(tmp)

            query = tmp / "query.fasta"

            output = tmp / "blast.tsv"

            SeqIO.write(
                SeqRecord(
                    Seq(sequence),
                    id="primer",
                    description="",
                ),
                query,
                "fasta",
            )

            self.runner.search(
                query,
                database,
                output,
            )

            return self.parser.parse(
                output,
            )

    def validate(
        self,
        sequence: str,
        database: str,
    ) -> float:

        hits = self.blast_hits(
            sequence,
            database,
        )

        return self.scorer.score(
            hits,
        )
