"""
Conserved region discovery.
"""

from pathlib import Path

from Bio import AlignIO

from primerforge.models.region import Region
from primerforge.alignment.coordinate_mapper import CoordinateMapper
from primerforge.analysis.primer_match import PrimerMatcher


class ConservedRegionFinder:
    """
    Find highly conserved windows in a multiple sequence alignment.
    """

    def __init__(self):

        self.last_statistics = {}

    def find(
        self,
        alignment: Path,
        window: int = 120,
        threshold: float = 95.0,
    ) -> list[Region]:

        aln = AlignIO.read(
            alignment,
            "fasta",
        )

        length = aln.get_alignment_length()

        regions = []

        for start in range(
            length - window + 1
        ):

            score = 0.0

            for column in range(
                start,
                start + window,
            ):

                bases = [
                    b
                    for b in aln[:, column]
                    if b != "-"
                ]

                if not bases:
                    continue

                most_common = max(
                    set(bases),
                    key=bases.count,
                )

                score += (
                    bases.count(most_common)
                    / len(bases)
                )

            score = (
                score / window
            ) * 100.0

            if score >= threshold:

                regions.append(
                    Region(
                        start=start,
                        end=start + window - 1,
                        score=round(
                            score,
                            2,
                        ),
                    )
                )

        return regions

    def best(
        self,
        alignment: Path,
        window: int = 120,
    ) -> Region | None:

        regions = self.find(
            alignment,
            window=window,
            threshold=0.0,
        )

        if not regions:
            return None

        return max(
            regions,
            key=lambda r: r.score,
        )

    def primer_conservation(
        self,
        primer: str,
        alignment: Path,
        start: int,
        max_mismatches: int = 1,
    ) -> float:
        """
        Estimate primer conservation using
        reference→alignment coordinate mapping
        and gap-aware primer matching.
        """

        aln = AlignIO.read(
            alignment,
            "fasta",
        )

        if len(aln) == 0:
            return 0.0

        #
        # Convert Primer3 reference coordinate
        # into alignment coordinate.
        #
        mapper = CoordinateMapper(
            alignment,
        )

        try:

            alignment_start = (
                mapper.reference_to_alignment(
                    start,
                )
            )

        except ValueError:

            #
            # Invalid coordinate.
            #
            self.last_statistics = {
                "matched_sequences": 0,
                "total_sequences": len(aln),
                "mean_identity": 0.0,
                "mean_mismatches": 0.0,
            }

            return 0.0

        matcher = PrimerMatcher(
            max_mismatches=max_mismatches,
        )

        matched = 0

        identities = []

        mismatches = []

        for record in aln:

            result = matcher.match(
                primer=primer,
                sequence=str(record.seq),
                start=alignment_start,
            )

            identities.append(
                result.identity,
            )

            mismatches.append(
                result.mismatches,
            )

            if result.matched:
                matched += 1

        conservation = (
            matched
            / len(aln)
            * 100.0
        )

        self.last_statistics = {

            "matched_sequences": matched,

            "total_sequences": len(aln),

            "mean_identity": round(
                sum(identities)
                / len(identities),
                2,
            ),

            "mean_mismatches": round(
                sum(mismatches)
                / len(mismatches),
                2,
            ),

            "conservation": round(
                conservation,
                2,
            ),
        }

        return round(
            conservation,
            2,
        )
