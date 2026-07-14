"""
Primer coverage analysis.
"""

from pathlib import Path

from Bio import SeqIO

from envoprimer.analysis.primer_match import PrimerMatcher
from envoprimer.coverage.window import CoordinateWindow
from envoprimer.models.pair import PrimerPair


class CoverageAnalyzer:
    """
    Estimate primer coverage across a sequence collection.

    Uses expected primer coordinates instead of scanning the
    entire alignment.
    """

    def __init__(
        self,
        max_mismatches: int = 0,
        tolerance: int = 3,
    ):

        self.matcher = PrimerMatcher(
            max_mismatches=max_mismatches,
        )

        self.window = CoordinateWindow(
            tolerance=tolerance,
        )

    def _matches(
        self,
        primer: str,
        sequence: str,
        expected_start: int,
    ) -> bool:

        for start in self.window.starts(
            expected=expected_start,
            sequence_length=len(sequence),
        ):

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

        forward_start = pair.forward.coordinate.start
        reverse_start = pair.reverse.coordinate.start

        for record in SeqIO.parse(
            fasta,
            "fasta",
        ):

            total += 1

            sequence = str(record.seq)

            forward = self._matches(
                pair.forward.sequence,
                sequence,
                forward_start,
            )

            reverse = self._matches(
                pair.reverse.sequence,
                sequence,
                reverse_start,
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
                forward_matches / total * 100,
                2,
            ),

            "reverse": round(
                reverse_matches / total * 100,
                2,
            ),

            "pair": round(
                pair_matches / total * 100,
                2,
            ),

            "forward_count": forward_matches,

            "reverse_count": reverse_matches,

            "pair_count": pair_matches,

            "total": total,

        }
