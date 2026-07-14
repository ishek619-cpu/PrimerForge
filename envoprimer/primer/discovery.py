"""
Complete primer discovery workflow.
"""

from pathlib import Path

from envoprimer.analysis.population import PopulationAnalyzer
from envoprimer.pipeline.validator import PrimerValidationPipeline
from envoprimer.primer.deduplicate import PrimerPairDeduplicator
from envoprimer.primer.pair import PrimerPairGenerator
from envoprimer.primer3.designer import Primer3Designer
from envoprimer.reference.mapper import CoordinateMapper
from envoprimer.reference.reference import Reference
from envoprimer.specificity.blastdb import BlastDatabase
from envoprimer.specificity.engine import SpecificityEngine


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
        region=None,
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
            template=reference.sequence,
            region=region,
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
        # Keep only the best candidates
        #
        pairs = pairs[: self.max_candidates]

        #
        # Population analysis
        #
        if mapper is not None:

            for pair in pairs:

                #
                # Convert coordinates ONLY for population analysis.
                # Do NOT overwrite primer coordinates.
                #
                forward_coordinate = mapper.reference_to_alignment(
                    pair.forward.coordinate,
                )

                reverse_coordinate = mapper.reference_to_alignment(
                    pair.reverse.coordinate,
                )

                forward_result = self.population.analyse(
                    primer=pair.forward.sequence,
                    alignment=reference.alignment,
                    coordinate=forward_coordinate,
                )

                reverse_result = self.population.analyse(
                    primer=pair.reverse.sequence,
                    alignment=reference.alignment,
                    coordinate=reverse_coordinate,
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
        # Validation
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
