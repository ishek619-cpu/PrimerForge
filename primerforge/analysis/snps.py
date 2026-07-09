"""
SNP discovery from multiple sequence alignments.
"""

from pathlib import Path

from Bio import AlignIO

from primerforge.models.snp import SNP


class SNPFinder:

    def find(
        self,
        alignment_file: Path,
        gap_threshold: float = 0.2,
    ) -> list[SNP]:

        alignment = AlignIO.read(
            alignment_file,
            "fasta",
        )

        nseq = len(alignment)

        length = alignment.get_alignment_length()

        snps = []

        for position in range(length):

            counts = {}

            gaps = 0

            for record in alignment:

                base = record.seq[position].upper()

                if base == "-":
                    gaps += 1
                    continue

                counts[base] = counts.get(base, 0) + 1

            if gaps / nseq > gap_threshold:
                continue

            if len(counts) <= 1:
                continue

            reference = max(
                counts,
                key=counts.get,
            )

            alternatives = [
                b
                for b in counts
                if b != reference
            ]

            snps.append(
                SNP(
                    position=position,
                    reference=reference,
                    alternatives=alternatives,
                    counts=counts,
                )
            )

        return snps

    def summary(
        self,
        alignment_file: Path,
    ):

        snps = self.find(
            alignment_file,
        )

        alignment = AlignIO.read(
            alignment_file,
            "fasta",
        )

        print()

        print(f"Sequences : {len(alignment)}")
        print(f"Alignment : {alignment.get_alignment_length()} bp")
        print(f"SNPs      : {len(snps)}")

        print()

        return snps
