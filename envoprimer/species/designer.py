"""
Species-specific primer designer.

Bridges the preprocessing workflow with the EnvoPrimer engine
and performs deduplication before expensive scoring.
"""

from __future__ import annotations

from Bio import SeqIO

from envoprimer.species.workflow import WorkflowResult

from envoprimer_v2.primer.engine import EnvoPrimerEngine

from envoprimer.conservation.engine import ConservationEngine
from envoprimer.specificity.engine import SpecificityEngine
from envoprimer.scoring.engine import ScoringEngine


class SpeciesPrimerDesigner:

    def __init__(self):

        self.engine = EnvoPrimerEngine()

        self.conservation = ConservationEngine()

        self.specificity = SpecificityEngine()

        self.scoring = ScoringEngine()

    ############################################################

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

    ############################################################

    def design(
        self,
        workflow: WorkflowResult,
    ):

        if not workflow.diagnostic_windows:

            raise RuntimeError(
                "No diagnostic windows were found."
            )

        print()
        print("=" * 60)
        print("ALIGNMENTS")
        print("=" * 60)

        target_alignment = self._load_alignment(
            workflow.alignment.target_alignment,
        )

        background_sequences = self._load_alignment(
            workflow.alignment.background_alignment,
        )

        print(f"Target sequences     : {len(target_alignment)}")
        print(f"Background sequences : {len(background_sequences)}")

        ############################################################
        # Generate primer pairs
        ############################################################

        pairs = []

        for window in workflow.diagnostic_windows:

            result = self.engine.design(

                window=window,

                diagnostic_snps=workflow.diagnostic_sites,

                background_sequences=background_sequences,

            )

            pairs.extend(result)

        ############################################################
        # REMOVE DUPLICATES
        ############################################################

        print()
        print("=" * 60)
        print("DEDUPLICATION")
        print("=" * 60)

        before = len(pairs)

        unique = {}

        for pair in pairs:

            key = (

                pair.forward.sequence,

                pair.reverse.sequence,

                pair.product_size,

            )

            #
            # Keep highest scoring duplicate
            #

            if key not in unique:

                unique[key] = pair

            else:

                old = unique[key]

                old_score = getattr(old, "primer3_penalty", 9999)

                new_score = getattr(pair, "primer3_penalty", 9999)

                if new_score < old_score:

                    unique[key] = pair

        pairs = list(unique.values())

        after = len(pairs)

        print(f"Before : {before}")
        print(f"After  : {after}")
        print(f"Removed: {before-after}")

        ############################################################
        # EARLY FILTER
        ############################################################

        print()
        print("=" * 60)
        print("EARLY FILTER")
        print("=" * 60)

        pairs.sort(

            key=lambda p: getattr(

                p,

                "primer3_penalty",

                9999,

            )

        )

        print(f"Before filtering : {len(pairs)}")

        MAX_PAIRS = 200

        pairs = pairs[:MAX_PAIRS]

        print(f"After filtering  : {len(pairs)}")

        ############################################################
        # BUILD BACKGROUND INDEX
        ############################################################

        self.specificity.build(
            background_sequences,
        )

        ############################################################
        # SCORING
        ############################################################

        print()
        print("=" * 60)
        print("SCORING")
        print("=" * 60)

        for pair in pairs:

            self.conservation.score(
                pair,
                target_alignment,
            )

            self.specificity.score(
                pair,
                background_sequences,
            )

            self.scoring.score(
                pair,
            )

        ############################################################

        pairs.sort(

            key=lambda p: p.final_score,

            reverse=True,

        )

        return pairs
