"""
PCR Engine.

Coordinates candidate searching, primer alignment,
binding evaluation and PCR product prediction.
"""

from __future__ import annotations

from envoprimer.specificity.index import BackgroundIndex
from envoprimer.specificity.search import CandidateSearcher

from envoprimer.pcr.alignment import PrimerAligner
from envoprimer.pcr.product import PCRProductPredictor
from envoprimer.pcr.evaluator import PCREvaluator


class PCREngine:
    """
    Exact in-silico PCR engine.
    """

    def __init__(
        self,
        seed_length: int = 18,
        minimum_product: int = 70,
        maximum_product: int = 200,
    ):

        self.index = BackgroundIndex(
            k=seed_length,
        )

        self.searcher = CandidateSearcher(
            seed_length=seed_length,
        )

        self.aligner = PrimerAligner()

        self.evaluator = PCREvaluator()

        self.predictor = PCRProductPredictor(
            minimum_product=minimum_product,
            maximum_product=maximum_product,
        )

    ############################################################

    def build(
        self,
        background_sequences: list[str],
    ):

        print()
        print("=" * 60)
        print("BUILDING PCR DATABASE")
        print("=" * 60)

        self.index.build(
            background_sequences,
        )

    ############################################################

    def _evaluate_hits(
        self,
        primer_sequence,
        hits,
    ):

        accepted = []

        for hit in hits:

            template = self.index.sequence(
                hit.sequence_id,
            )

            ####################################################
            # CandidateSearcher indexes the LAST seed_length
            # bases of the primer.
            #
            # Convert the seed coordinate back to the start
            # of the FULL primer.
            ####################################################

            primer_start = hit.position - (
                len(primer_sequence)
                - self.searcher.seed_length
            )

            binding = self.aligner.align(

                primer=primer_sequence,

                template=template,

                sequence_id=hit.sequence_id,

                start=primer_start,

                strand=hit.strand,

            )

            ####################################################
            # DEBUG
            ####################################################

            print(
                f"Seq={binding.sequence_id} "
                f"strand={binding.strand} "
                f"start={binding.start} "
                f"end={binding.end} "
                f"identity={binding.identity:.2f} "
                f"mismatches={binding.mismatches} "
                f"terminal={binding.terminal_mismatches} "
                f"score={binding.score:.2f}"
            )

            passed = self.evaluator.passes(
                binding,
            )

            print(
                "PASS"
                if passed
                else "FAIL"
            )

            if passed:

                accepted.append(
                    binding,
                )

        return accepted

    ############################################################

    def predict(
        self,
        pair,
    ):

        ########################################################
        # Candidate search
        ########################################################

        forward_hits = self.searcher.search(

            pair.forward.sequence,

            self.index,

        )

        reverse_hits = self.searcher.search(

            pair.reverse.sequence,

            self.index,

        )

        print()
        print("=" * 60)
        print("CANDIDATES")
        print("=" * 60)

        print(
            f"Forward candidates : {len(forward_hits)}"
        )

        print(
            f"Reverse candidates : {len(reverse_hits)}"
        )

        ########################################################
        # Alignment
        ########################################################

        print()
        print("=" * 60)
        print("FORWARD ALIGNMENTS")
        print("=" * 60)

        forward_bindings = self._evaluate_hits(

            pair.forward.sequence,

            forward_hits,

        )

        print()
        print("=" * 60)
        print("REVERSE ALIGNMENTS")
        print("=" * 60)

        reverse_bindings = self._evaluate_hits(

            pair.reverse.sequence,

            reverse_hits,

        )

        print()
        print("=" * 60)
        print("ACCEPTED BINDINGS")
        print("=" * 60)

        print(
            f"Forward : {len(forward_bindings)}"
        )

        print(
            f"Reverse : {len(reverse_bindings)}"
        )

        ########################################################
        # PCR prediction
        ########################################################

        products = self.predictor.predict(

            forward_bindings,

            reverse_bindings,

        )

        print()

        print(
            f"PCR Products : {len(products)}"
        )

        return products
