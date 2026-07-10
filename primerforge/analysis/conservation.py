"""
Conserved region discovery.
"""

from pathlib import Path

from Bio import AlignIO

from primerforge.models.region import Region


class ConservedRegionFinder:
    """
    Find highly conserved windows in a multiple sequence alignment.
    """

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
