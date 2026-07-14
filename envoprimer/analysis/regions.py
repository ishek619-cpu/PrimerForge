"""
Candidate region discovery.
"""

from envoprimer.models.conservation import ConservationProfile
from envoprimer.models.region import Region
from envoprimer.models.snp import SNP


class RegionFinder:
    """
    Discover candidate primer regions around informative SNPs.
    """

    def __init__(
        self,
        flank: int = 125,
        minimum_average_conservation: float = 0.70,
    ):

        self.flank = flank
        self.minimum_average_conservation = minimum_average_conservation

    def find(
        self,
        profile: ConservationProfile,
        snps: list[SNP],
    ) -> list[Region]:

        regions = []

        for snp in snps:

            start = max(
                0,
                snp.position - self.flank,
            )

            end = min(
                profile.alignment_length - 1,
                snp.position + self.flank,
            )

            scores = profile.scores[start:end + 1]

            if not scores:
                continue

            average = sum(scores) / len(scores)

            if average < self.minimum_average_conservation:
                continue

            regions.append(
                Region(
                    snp_position=snp.position,
                    start=start,
                    end=end,
                    length=end - start + 1,
                )
            )

        return regions

    def summary(
        self,
        profile: ConservationProfile,
        snps: list[SNP],
    ):

        regions = self.find(
            profile,
            snps,
        )

        print()

        print(f"SNPs              : {len(snps)}")
        print(f"Candidate regions : {len(regions)}")

        print()

        return regions
