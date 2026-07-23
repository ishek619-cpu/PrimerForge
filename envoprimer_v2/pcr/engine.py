"""
EnvoPrimer v2 PCR Engine.
"""

from __future__ import annotations

from envoprimer_v2.pcr.search import SeedSearcher
from envoprimer_v2.pcr.window import WindowExtractor
from envoprimer_v2.pcr.local_aligner import LocalAligner
from envoprimer_v2.pcr.evaluator import BindingEvaluator
from envoprimer_v2.pcr.product import ProductPredictor


class PCREngine:
    """
    End-to-end in silico PCR engine.
    """

    def __init__(
        self,
        seed_length: int = 18,
        flank: int = 25,
    ):

        self.searcher = SeedSearcher(seed_length)
        self.window = WindowExtractor(flank)
        self.aligner = LocalAligner()
        self.evaluator = BindingEvaluator()
        self.predictor = ProductPredictor()

    def build(self, sequences: list[str]):

        self.sequences = [s.upper() for s in sequences]
        self.searcher.build(self.sequences)

    def _find_sites(
        self,
        primer: str,
        strand: str,
    ):

        accepted = []
        seen = set()

        # Searcher now handles reverse-complement logic.
        candidates = self.searcher.search(
            primer,
            reverse=(strand == "-"),
        )

        for sequence_id, seed_start in candidates:

            sequence = self.sequences[sequence_id]

            window, window_start = self.window.extract(
                sequence=sequence,
                seed_start=seed_start,
                seed_length=self.searcher.seed_length,
            )

            site = self.aligner.align(
                primer=primer,
                window=window,
                sequence_id=sequence_id,
                strand=strand,
                window_start=window_start,
            )

            if site is None:
                continue

            if not self.evaluator.accept(site):
                continue

            key = (
                site.sequence_id,
                site.strand,
                site.start,
                site.end,
            )

            if key in seen:
                continue

            seen.add(key)
            accepted.append(site)

        return accepted

    def predict(
        self,
        forward_primer: str,
        reverse_primer: str,
    ):

        forward_sites = self._find_sites(
            forward_primer,
            "+",
        )

        reverse_sites = self._find_sites(
            reverse_primer,
            "-",
        )

        return self.predictor.predict(
            forward_sites,
            reverse_sites,
        )
