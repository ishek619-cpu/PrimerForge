"""
Primer discovery engine.

Designs, pairs, deduplicates and ranks primers.
"""

from pathlib import Path

from Bio import SeqIO

from primerforge.primer3.designer import Primer3Designer
from primerforge.primer.pair import PrimerPairGenerator
from primerforge.primer.pairscore import PrimerPairScorer
from primerforge.primer.deduplicate import PrimerPairDeduplicator


class PrimerDiscovery:
    """
    Complete primer discovery workflow.
    """

    def __init__(self):

        self.designer = Primer3Designer()

        self.generator = PrimerPairGenerator(
            min_product=80,
            max_product=250,
        )

        self.scorer = PrimerPairScorer()

        self.deduplicator = PrimerPairDeduplicator()

    def discover(
        self,
        fasta: Path,
    ):

        record = next(
            SeqIO.parse(
                fasta,
                "fasta",
            )
        )

        forward, reverse = self.designer.design(
            str(record.seq)
        )

        pairs = self.generator.generate(
            forward,
            reverse,
        )

        pairs = self.deduplicator.deduplicate(
            pairs,
        )

        ranked = self.scorer.rank(
            pairs,
        )

        return ranked
