"""
EnvoPrimer Core Pipeline.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from envoprimer.cache.engine import CacheEngine

from envoprimer.species.workflow import SpeciesWorkflow
from envoprimer.species.designer import SpeciesPrimerDesigner

from envoprimer.conservation.engine import ConservationEngine
from envoprimer.specificity.engine import SpecificityEngine
from envoprimer.scoring.engine import ScoringEngine


class EnvoPrimerPipeline:
    """
    Main EnvoPrimer pipeline.
    """

    def __init__(
        self,
        email: str,
        api_key: str | None = None,
    ):

        self.workflow = SpeciesWorkflow(
            email=email,
            api_key=api_key,
        )

        self.designer = SpeciesPrimerDesigner()

        self.conservation = ConservationEngine()

        self.specificity = SpecificityEngine()

        self.scoring = ScoringEngine()

        self.cache = CacheEngine()

    def _load_alignment(
        self,
        fasta,
    ) -> list[str]:

        return [

            str(record.seq).upper()

            for record in SeqIO.parse(
                fasta,
                "fasta",
            )

        ]

    def run(
        self,
        species: str,
        marker: str,
        output_directory: str | Path,
    ):

        output_directory = Path(output_directory)

        print()
        print("=" * 70)
        print("ENVOPRIMER")
        print("=" * 70)

        ##########################################################
        # CACHE
        ##########################################################

        if self.cache.exists(
            species,
            marker,
        ):

            print()
            print("=" * 70)
            print("CACHE FOUND")
            print("=" * 70)

            metadata = self.cache.load_metadata(
                species,
                marker,
            )

            if metadata is not None:

                print(
                    f"Cached Primer Pairs : {metadata.get('primer_pairs', 0)}"
                )

        else:

            print()
            print("=" * 70)
            print("CACHE NOT FOUND")
            print("=" * 70)

        ##########################################################
        # WORKFLOW
        ##########################################################

        workflow = self.workflow.run(

            species=species,

            marker=marker,

            output_directory=output_directory,

        )

        ##########################################################
        # ALIGNMENT
        ##########################################################

        target_alignment = self._load_alignment(

            workflow.alignment.target_alignment

        )

        ##########################################################
        # PRIMER DESIGN
        ##########################################################

        pairs = self.designer.design(

            workflow,

        )

        ##########################################################
        # SCORING
        ##########################################################

        for pair in pairs:

            #
            # Population conservation
            #

            self.conservation.score(

                pair,

                target_alignment,

            )

            #
            # Species specificity
            #

            self.specificity.score(

                pair,

                [],

            )

            #
            # Final score
            #

            self.scoring.score(

                pair,

            )

        ##########################################################
        # SORT
        ##########################################################

        pairs.sort(

            key=lambda p: p.final_score,

            reverse=True,

        )

        ##########################################################
        # SAVE CACHE
        ##########################################################

        self.cache.save_metadata(

            species,

            marker,

            {

                "species": species,

                "marker": marker,

                "primer_pairs": len(pairs),

            },

        )

        ##########################################################
        # FINISHED
        ##########################################################

        print()
        print("=" * 70)
        print("PIPELINE COMPLETE")
        print("=" * 70)

        print(
            f"Primer pairs : {len(pairs)}"
        )

        return pairs
# ----------------------------------------------------------------------
# Backward compatibility
# ----------------------------------------------------------------------

Pipeline = EnvoPrimerPipeline
