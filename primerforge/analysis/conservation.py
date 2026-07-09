"""
Alignment conservation analysis.
"""

from pathlib import Path

from Bio import AlignIO

from primerforge.models.conservation import ConservationProfile


class ConservationAnalyzer:
    """
    Analyze sequence conservation in a multiple sequence alignment.
    """

    def calculate(
        self,
        alignment_file: Path,
    ) -> ConservationProfile:

        alignment = AlignIO.read(
            alignment_file,
            "fasta",
        )

        nseq = len(alignment)

        length = alignment.get_alignment_length()

        scores = []

        for position in range(length):

            counts = {}

            for record in alignment:

                base = record.seq[position].upper()

                if base == "-":
                    continue

                counts[base] = counts.get(base, 0) + 1

            if not counts:
                scores.append(0.0)
                continue

            scores.append(
                max(counts.values()) / nseq
            )

        return ConservationProfile(
            alignment_length=length,
            number_of_sequences=nseq,
            scores=scores,
        )

    def summary(
        self,
        alignment_file: Path,
    ) -> ConservationProfile:

        profile = self.calculate(
            alignment_file,
        )

        print()

        print(
            f"Sequences : {profile.number_of_sequences}"
        )

        print(
            f"Length    : {profile.alignment_length}"
        )

        print(
            f"Average conservation : {profile.average:.3f}"
        )

        print()

        return profile
