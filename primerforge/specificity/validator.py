"""
Primer specificity validation.
"""

from pathlib import Path
from tempfile import TemporaryDirectory

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from primerforge.specificity.blast import BlastRunner
from primerforge.specificity.parser import BlastParser
from primerforge.specificity.scorer import SpecificityScorer


class SpecificityValidator:

    def __init__(self):

        self.runner = BlastRunner()

        self.parser = BlastParser()

        self.scorer = SpecificityScorer()

    def validate(
        self,
        sequence: str,
        database: str,
    ) -> float:

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

            hits = self.parser.parse(
                output,
            )

            return self.scorer.score(
                hits,
            )
