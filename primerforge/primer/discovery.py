"""
Complete primer discovery workflow.
"""

from pathlib import Path

from primerforge.analysis.population import PopulationAnalyzer
from primerforge.pipeline.validator import PrimerValidationPipeline
from primerforge.primer.deduplicate import PrimerPairDeduplicator
from primerforge.primer.pair import PrimerPairGenerator
from primerforge.primer3.designer import Primer3Designer
from primerforge.reference.mapper import CoordinateMapper
from primerforge.reference.reference import Reference
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
        max_candidates: int = 20,
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

        self.specificity = None

        #
        # Only validate the best N candidates.
        #
        self.max_candidates = max_candidates

    def discover(
        self,
        reference_fasta: Path,
        specificity_engine: SpecificityEngine | None = None,
        alignment_fasta: Path | None = None,
    ):

        reference = Reference(
            fasta=reference_fasta,
            alignment=alignment_fasta,
        )

        mapper = None

        if reference.has_alignment:

            mapper = CoordinateMapper(
                reference,
            )

        #
        # Primer3
        #
        forward, reverse = self.designer.design(
            reference.sequence,
        )

        #
        # Generate pairs
        #
        pairs = self.generator.generate(
            forward,
            reverse,
        )

        #
        # Remove duplicates
        #
        pairs = self.deduplicator.deduplicate(
            pairs,
        )

        #
        # Generator already sorts pairs by quality.
        # Keep only the top candidates.
        #
        pairs = pairs[: self.max_candidates]

        #
        # Population analysis
        #
        if mapper is not None:

            for pair in pairs:

                forward_coordinate = mapper.reference_to_alignment(
                    pair.forward.coordinate,
                )

                reverse_coordinate = mapper.reference_to_alignment(
                    pair.reverse.coordinate,
                )

                forward_result = self.population.analyse(
                    primer=pair.forward.sequence,
                    alignment=reference.alignment,
                    start=forward_coordinate.start,
                )

                reverse_result = self.population.analyse(
                    primer=pair.reverse.sequence,
                    alignment=reference.alignment,
                    start=reverse_coordinate.start,
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
        # Expensive validation only for the top candidates.
        #
        engine = specificity_engine or self.specificity

        pairs = self.validator.validate(
            pairs=pairs,
            reference_fasta=reference.fasta,
            alignment_fasta=reference.alignment,
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
