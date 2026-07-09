"""
Degenerate primer generator using IUPAC ambiguity codes.
"""

from Bio import AlignIO


IUPAC = {
    frozenset({"A"}): "A",
    frozenset({"C"}): "C",
    frozenset({"G"}): "G",
    frozenset({"T"}): "T",
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "G"}): "V",
    frozenset({"A", "C", "G", "T"}): "N",
}


class DegeneratePrimerGenerator:
    """
    Convert an alignment region into a degenerate primer.
    """

    def generate(
        self,
        alignment_file,
        start,
        length,
    ):

        alignment = AlignIO.read(
            alignment_file,
            "fasta",
        )

        primer = []

        for position in range(
            start,
            start + length,
        ):

            bases = set()

            for record in alignment:

                base = record.seq[position].upper()

                if base != "-":

                    bases.add(base)

            primer.append(
                IUPAC[
                    frozenset(bases)
                ]
            )

        return "".join(primer)
