"""
Diagnostic window builder.

Creates candidate Primer3 windows around clusters of
diagnostic SNPs.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer_v2.diagnostics.snp import DiagnosticSNP


@dataclass(slots=True)
class DiagnosticWindow:

    start: int

    end: int

    sequence: str

    diagnostic_snps: list[DiagnosticSNP]

    conservation: float = 0.0

    score: float = 0.0

    @property
    def length(self):

        return self.end - self.start + 1


class DiagnosticWindowBuilder:

    def __init__(

        self,

        minimum_size: int = 120,

        maximum_size: int = 200,

    ):

        self.minimum_size = minimum_size

        self.maximum_size = maximum_size

    def build(

        self,

        reference_sequence: str,

        snps: list[DiagnosticSNP],

    ) -> list[DiagnosticWindow]:

        windows = []

        if not snps:

            return windows

        used = set()

        for snp in snps:

            start = max(

                0,

                snp.alignment_position

                - self.minimum_size // 2,

            )

            end = min(

                len(reference_sequence),

                start + self.minimum_size,

            )

            start = max(

                0,

                end - self.minimum_size,

            )

            key = (

                start,

                end,

            )

            if key in used:

                continue

            used.add(key)

            contained = [

                s

                for s in snps

                if start <= s.alignment_position <= end

            ]

            windows.append(

                DiagnosticWindow(

                    start=start,

                    end=end,

                    sequence=reference_sequence[start:end],

                    diagnostic_snps=contained,

                    score=float(len(contained)),

                )

            )

        windows.sort(

            key=lambda w: (

                -len(w.diagnostic_snps),

                w.start,

            )

        )

        return windows
