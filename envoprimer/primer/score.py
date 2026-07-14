"""
Primer scoring utilities.
"""

from envoprimer.models.primer import Primer


class PrimerScorer:

    def gc_percent(
        self,
        sequence: str,
    ) -> float:

        sequence = sequence.upper()

        gc = sequence.count("G") + sequence.count("C")

        return gc * 100 / len(sequence)

    def wallace_tm(
        self,
        sequence: str,
    ) -> float:

        sequence = sequence.upper()

        a = sequence.count("A")
        t = sequence.count("T")
        g = sequence.count("G")
        c = sequence.count("C")

        return 2 * (a + t) + 4 * (g + c)

    def has_gc_clamp(
        self,
        sequence: str,
    ) -> bool:

        return sequence[-1] in "GC"

    def longest_homopolymer(
        self,
        sequence: str,
    ) -> int:

        longest = 1
        current = 1

        for i in range(1, len(sequence)):

            if sequence[i] == sequence[i - 1]:

                current += 1

                longest = max(
                    longest,
                    current,
                )

            else:

                current = 1

        return longest

    def evaluate(
        self,
        primer: Primer,
    ) -> Primer:

        primer.gc = self.gc_percent(
            primer.sequence,
        )

        primer.tm = self.wallace_tm(
            primer.sequence,
        )

        primer.gc_clamp = self.has_gc_clamp(
            primer.sequence,
        )

        primer.homopolymer = (
            self.longest_homopolymer(
                primer.sequence,
            ) >= 5
        )

        return primer
