"""
EnvoPrimer complete design pipeline.
"""

from __future__ import annotations

from envoprimer_v2.models.alignment import Alignment
from envoprimer_v2.diagnostics.snp import DiagnosticSNPFinder
from envoprimer_v2.diagnostics.conservation import ConservationAnalyzer
from envoprimer_v2.windows.builder import DiagnosticWindowBuilder
from envoprimer_v2.primer.engine import EnvoPrimerEngine


class DesignPipeline:

    def __init__(self):

        self.snp = DiagnosticSNPFinder()

        self.conservation = ConservationAnalyzer()

        self.windows = DiagnosticWindowBuilder()

        self.engine = EnvoPrimerEngine()

    def run(

        self,

        reference_sequence: str,

        target_alignment_fasta,

        background_alignment_fasta,

        background_sequences,

    ):

        print()

        print("=" * 70)

        print("Building alignments")

        print("=" * 70)

        target = Alignment.from_fasta(

            target_alignment_fasta,

        )

        background = Alignment.from_fasta(

            background_alignment_fasta,

        )

        print(

            f"Target sequences     : {target.count}"

        )

        print(

            f"Background sequences : {background.count}"

        )

        print()

        print("=" * 70)

        print("Finding diagnostic SNPs")

        print("=" * 70)

        diagnostic_snps = self.snp.find(

            target,

            background,

        )

        print(

            f"Diagnostic SNPs : {len(diagnostic_snps)}"

        )

        print()

        print("=" * 70)

        print("Building diagnostic windows")

        print("=" * 70)

        windows = self.windows.build(

            reference_sequence,

            diagnostic_snps,

        )

        print(

            f"Candidate windows : {len(windows)}"

        )

        print()

        conserved = self.conservation.analyse(

            target,

        )

        all_pairs = []

        for window in windows:

            window.conservation = self.conservation.window_score(

                conserved,

                window.start,

                window.end,

            )

            window.score += window.conservation

            pairs = self.engine.design(

                window=window,

                diagnostic_snps=diagnostic_snps,

                background_sequences=background_sequences,

            )

            all_pairs.extend(

                pairs,

            )

        all_pairs.sort(

            key=lambda pair: pair.final_score,

            reverse=True,

        )

        return all_pairs
