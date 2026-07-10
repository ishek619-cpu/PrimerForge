"""
Complete primer discovery workflow.
"""

from pathlib import Path

from Bio import SeqIO

from primerforge.analysis.population import PopulationAnalyzer
from primerforge.pipeline.validator import PrimerValidationPipeline
from primerforge.primer.deduplicate import PrimerPairDeduplicator
from primerforge.primer.pair import PrimerPairGenerator
from primerforge.primer3.designer import Primer3Designer
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

        self.population = PopulationAnalyzer()

        self.validator = PrimerValidationPipeline()

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

        #
        # Load reference sequence
        #
        record = next(
            SeqIO.parse(
                reference_fasta,
                "fasta",
            )
        )

        sequence = str(record.seq)

        #
        # Design primers
        #
        forward, reverse = self.designer.design(
            sequence,
        )

        #
        # Generate candidate primer pairs
        #
        pairs = self.generator.generate(
            forward,
            reverse,
        )

        #
        # Remove duplicate pairs
        #
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
        # Complete validation pipeline
        #
        engine = specificity_engine or self.specificity

        pairs = self.validator.validate(
            pairs=pairs,
            reference_fasta=reference_fasta,
            alignment_fasta=alignment_fasta,
            specificity_engine=engine,
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
            reference_fasta=reference_fasta,
            alignment_fasta=alignment_fasta,
        )
