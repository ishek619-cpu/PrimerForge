"""
Complete primer discovery workflow.
"""

from pathlib import Path

from Bio import SeqIO

from primerforge.primer3.designer import Primer3Designer
from primerforge.primer.pair import PrimerPairGenerator
from primerforge.primer.pairscore import PrimerPairScorer
from primerforge.primer.deduplicate import PrimerPairDeduplicator

from primerforge.analysis.population import PopulationAnalyzer

from primerforge.specificity.blastdb import BlastDatabase
from primerforge.specificity.engine import SpecificityEngine


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

        self.population = PopulationAnalyzer()

        self.blastdb = BlastDatabase()

        #
        # Optional specificity engine.
        #
        self.specificity = None

    def discover(
        self,
        reference_fasta: Path,
        specificity_engine: SpecificityEngine | None = None,
        alignment_fasta: Path | None = None,
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
        # Population analysis
        #
        if alignment_fasta is not None:

            for pair in pairs:

                forward_result = self.population.analyse(
                    primer=pair.forward.sequence,
                    alignment=alignment_fasta,
                    start=pair.forward.start,
                )

                reverse_result = self.population.analyse(
                    primer=pair.reverse.sequence,
                    alignment=alignment_fasta,
                    start=pair.reverse.start,
                )

                pair.population_conservation = min(
                    forward_result["population_conservation"],
                    reverse_result["population_conservation"],
                )

                pair.population_coverage = min(
                    forward_result["population_coverage"],
                    reverse_result["population_coverage"],
                )

                pair.population_result = {
                    "forward": forward_result,
                    "reverse": reverse_result,
                }

        #
        # Species specificity
        #
        engine = specificity_engine or self.specificity

        if engine is not None:

            blast_dir = Path(
                "results/blast"
            )

            for pair in pairs:

                engine.evaluate_pair(
                    pair,
                    blast_dir,
                )

        #
        # Initial ranking
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

        database_prefix = str(
            database_prefix,
        )

        self.blastdb.build(
            alignment_fasta,
            database_prefix,
        )

        return self.discover(
            reference_fasta,
            alignment_fasta=alignment_fasta,
        )
