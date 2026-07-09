"""
PrimerForge Diagnostic Region Discovery Engine
"""

from pathlib import Path

from Bio import AlignIO


class RegionFinder:
    """
    Finds conserved regions surrounding species-specific SNPs.
    """

    def __init__(self):
        pass

    def find_candidate_regions(
        self,
        alignment_file: Path,
        snps: list,
        flank: int = 25,
    ):

        alignment = AlignIO.read(
            alignment_file,
            "fasta",
        )

        length = alignment.get_alignment_length()

        candidates = []

        for snp in snps:

            position = snp["position"] - 1

            start = max(0, position - flank)

            end = min(length, position + flank + 1)

            candidates.append(
                {
                    "snp_position": position + 1,
                    "start": start + 1,
                    "end": end,
                    "length": end - start,
                }
            )

        return candidates

    def summary(
        self,
        alignment_file: Path,
        snps: list,
    ):

        candidates = self.find_candidate_regions(
            alignment_file,
            snps,
        )

        print()

        print(f"Candidate regions : {len(candidates)}")

        print()

        for region in candidates[:20]:

            print(region)
