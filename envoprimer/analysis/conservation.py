"""
Conserved region discovery.
"""

from pathlib import Path

from Bio import AlignIO

from envoprimer.analysis.primer_match import PrimerMatcher
from envoprimer.models.region import Region


class ConservedRegionFinder:

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

        for start in range(length - window + 1):

            score = 0.0

            for column in range(start, start + window):

                bases = [
                    b
                    for b in aln[:, column]
                    if b != "-"
                ]

                if not bases:
                    continue

                common = max(
                    set(bases),
                    key=bases.count,
                )

                score += (
                    bases.count(common)
                    / len(bases)
                )

            score = score / window * 100

            if score >= threshold:

                regions.append(

                    Region(
                        start=start,
                        end=start + window - 1,
                        score=round(score, 2),
                    )

                )

        return regions

    def best(
        self,
        alignment: Path,
        window: int = 120,
    ):

        regions = self.find(
            alignment,
            window,
            threshold=0,
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

        aln = AlignIO.read(
            alignment,
            "fasta",
        )

        matcher = PrimerMatcher(
            max_mismatches=max_mismatches,
        )

        matched = 0

        identities = []
        mismatches = []

        printed = False

        for i, record in enumerate(aln):

            result = matcher.match(
                primer=primer,
                sequence=str(record.seq),
                start=start,
            )

            if not printed:

                print("\n================ TRACE ================")
                print("Sequence :", record.id)
                print("Primer   :", primer)
                print("Start    :", start)
                print("Coord    :", result.coordinate)
                print("Extract  :", result.aligned_sequence)
                print("Identity :", result.identity)
                print("Mismatch :", result.mismatches)
                print("Matched? :", result.matched)
                print("=======================================\n")

                printed = True

            identities.append(result.identity)
            mismatches.append(result.mismatches)

            if result.matched:
                matched += 1

        conservation = (
            matched
            / len(aln)
            * 100
        )

        self.last_statistics = {

            "matched_sequences": matched,
            "total_sequences": len(aln),
            "mean_identity": round(
                sum(identities) / len(identities),
                2,
            ),
            "mean_mismatches": round(
                sum(mismatches) / len(mismatches),
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
