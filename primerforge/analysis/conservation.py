"""
Alignment conservation analysis.
"""

from pathlib import Path

from Bio import AlignIO


class ConservationAnalyzer:
    """
    Calculate conservation across a multiple sequence alignment.
    """

    def load_alignment(self, alignment_file: Path):

        return AlignIO.read(
            alignment_file,
            "fasta",
        )

    def conservation_scores(
        self,
        alignment,
    ) -> list[float]:

        scores = []

        length = alignment.get_alignment_length()

        nseq = len(alignment)

        for column in range(length):

            counts = {}

            for record in alignment:

                base = record.seq[column]

                if base == "-":
                    continue

                counts[base] = counts.get(base, 0) + 1

            if not counts:

                scores.append(0.0)

                continue

            maximum = max(counts.values())

            scores.append(maximum / nseq)

        return scores

    def summary(
        self,
        alignment_file: Path,
    ):

        alignment = self.load_alignment(
            alignment_file
        )

        scores = self.conservation_scores(
            alignment
        )

        print()

        print(f"Sequences : {len(alignment)}")
        print(f"Length    : {alignment.get_alignment_length()}")

        print(
            f"Average conservation : {sum(scores)/len(scores):.3f}"
        )

        print()
