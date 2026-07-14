"""
Primer discovery from conserved candidate regions.
"""

from pathlib import Path

from Bio import SeqIO

from envoprimer.analysis.conservation import ConservationAnalyzer
from envoprimer.analysis.regions import RegionFinder
from envoprimer.analysis.snps import SNPFinder

from envoprimer.primer3.designer import Primer3Designer
from envoprimer.primer.deduplicate import PrimerPairDeduplicator
from envoprimer.primer.pair import PrimerPairGenerator
from envoprimer.primer.pairscore import PrimerPairScorer


class RegionPrimerDiscovery:

    def __init__(self):

        self.conservation = ConservationAnalyzer()

        self.snpfinder = SNPFinder()

        self.regionfinder = RegionFinder()

        self.designer = Primer3Designer()

        self.generator = PrimerPairGenerator(
            min_product=80,
            max_product=250,
        )

        self.scorer = PrimerPairScorer()

        self.deduplicator = PrimerPairDeduplicator()

    def discover(
        self,
        alignment_file: Path,
        reference_fasta: Path,
    ):

        profile = self.conservation.calculate(
            alignment_file,
        )

        snps = self.snpfinder.find(
            alignment_file,
        )

        regions = self.regionfinder.find(
            profile,
            snps,
        )

        reference = next(
            SeqIO.parse(
                reference_fasta,
                "fasta",
            )
        )

        sequence = str(reference.seq)

        all_pairs = []

        for region in regions:

            start = max(0, region.start)

            end = min(
                len(sequence),
                region.end,
            )

            template = sequence[start:end]

            if len(template) < 80:
                continue

            try:

                forward, reverse = self.designer.design(
                    template
                )

            except Exception:
                continue

            pairs = self.generator.generate(
                forward,
                reverse,
            )

            all_pairs.extend(
                pairs,
            )

        all_pairs = self.deduplicator.deduplicate(
            all_pairs,
        )

        ranked = self.scorer.rank(
            all_pairs,
        )

        return ranked
