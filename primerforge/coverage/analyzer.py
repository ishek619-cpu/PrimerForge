"""
Primer coverage analysis.
"""

from pathlib import Path

from Bio import SeqIO

from primerforge.analysis.primer_match import PrimerMatcher
from primerforge.models.pair import PrimerPair


class CoverageAnalyzer:
    """
    Estimate primer coverage across a sequence collection.
    """

    def __init__(
        self,
        max_mismatches: int = 0,
    ):

        self.matcher = PrimerMatcher(
            max_mismatches=max_mismatches,
        )

    def _matches(
        self,
        primer: str,
        sequence: str,
    ) -> bool:

        #
        # Search every possible alignment position.
        #
        for start in range(len(sequence)):

            result = self.matcher.match(
                primer=primer,
                sequence=sequence,
                start=start,
            )

            if result.matched:

                return True

        return False

    def analyse(
        self,
        pair: PrimerPair,
        fasta: Path,
    ) -> dict:

        total = 0

        forward_matches = 0

        reverse_matches = 0

        pair_matches = 0

        for record in SeqIO.parse(
            fasta,
            "fasta",
        ):

            total += 1

            sequence = str(
                record.seq,
            )

            forward = self._matches(
                pair.forward.sequence,
                sequence,
            )

            reverse = self._matches(
                pair.reverse.sequence,
                sequence,
            )

            if forward:

                forward_matches += 1

            if reverse:

                reverse_matches += 1

            if forward and reverse:

                pair_matches += 1

        if total == 0:

            return {

                "forward": 0.0,

                "reverse": 0.0,

                "pair": 0.0,

                "forward_count": 0,

                "reverse_count": 0,

                "pair_count": 0,

                "total": 0,

            }

        return {

            "forward": round(
                forward_matches
                / total
                * 100,
                2,
            ),

            "reverse": round(
                reverse_matches
                / total
                * 100,
                2,
            ),

            "pair": round(
                pair_matches
                / total
                * 100,
                2,
            ),

            "forward_count": forward_matches,

            "reverse_count": reverse_matches,

            "pair_count": pair_matches,

            "total": total,

        }
