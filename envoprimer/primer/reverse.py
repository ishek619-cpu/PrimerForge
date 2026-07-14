"""
Reverse complement utilities.
"""

from envoprimer.models.primer import Primer


class ReverseComplement:

    _TABLE = str.maketrans(
        "ATCGNatcgn",
        "TAGCNtagcn",
    )

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """
        Return the reverse complement of a DNA sequence.
        """

        return sequence.translate(
            ReverseComplement._TABLE
        )[::-1]

    def reverse_primer(
        self,
        primer: Primer,
    ) -> Primer:
        """
        Convert a forward primer into a reverse primer.
        """

        return Primer(
            sequence=self.reverse_complement(
                primer.sequence
            ),
            start=primer.start,
            end=primer.end,
            strand="-",
            length=primer.length,
            tm=primer.tm,
            gc=primer.gc,
            gc_clamp=primer.gc_clamp,
            homopolymer=primer.homopolymer,
            hairpin_score=primer.hairpin_score,
            self_dimer_score=primer.self_dimer_score,
            score=primer.score,
            notes=list(primer.notes),
        )
