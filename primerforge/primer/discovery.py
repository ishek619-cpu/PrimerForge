"""
Complete primer discovery workflow.
"""

from pathlib import Path

from Bio import SeqIO

from primerforge.primer3.designer import Primer3Designer
from primerforge.primer.pair import PrimerPairGenerator
from primerforge.primer.pairscore import PrimerPairScorer
from primerforge.primer.deduplicate import PrimerPairDeduplicator
from primerforge.specificity.blastdb import BlastDatabase


class PrimerDiscovery:
    """
    Complete primer discovery workflow.
    """

    def __init__(
        self,
        min_product: int = 80,
        max_product: int = 250,
    ):

        self.designer = Primer3Designer()

        self.generator = PrimerPairGenerator(
            min_product=min_product,
            max_product=max_product,
        )

        self.deduplicator = PrimerPairDeduplicator()

        self.scorer = PrimerPairScorer()

        self.blastdb = BlastDatabase()

    def discover(
        self,
        reference_fasta: Path,
    ):

        record = next(
            SeqIO.parse(
                reference_fasta,
                "fasta",
            )
        )

        sequence = str(record.seq)

        forward, reverse = self.designer.design(
            sequence,
        )

        pairs = self.generator.generate(
            forward,
            reverse,
        )

        pairs = self.deduplicator.deduplicate(
            pairs,
        )

        #
        # Initial ranking only.
        # Validation later computes the real score.
        #
        for pair in pairs:

            validation = {
                "thermo": 100.0,
                "conservation": 100.0,
                "snp": 100.0,
            }

            self.scorer.score(
                pair,
                validation,
            )

        pairs = self.scorer.rank(
            pairs,
        )

        return pairs

    def discover_from_database(
        self,
        reference_fasta: Path,
        alignment_fasta: Path,
        database_prefix: str = "data/blast/primerforge",
    ):

        database_prefix = str(database_prefix)

        self.blastdb.build(
            alignment_fasta,
            database_prefix,
        )

        return self.discover(
            reference_fasta,
        )
