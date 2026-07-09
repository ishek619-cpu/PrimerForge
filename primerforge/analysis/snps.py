"""
PrimerForge Species-Specific SNP Discovery Engine
"""

from pathlib import Path

from Bio import AlignIO


class SNPFinder:

    def __init__(self):
        pass

    def find_variable_sites(
        self,
        alignment_file: Path,
    ):

        alignment = AlignIO.read(
            alignment_file,
            "fasta",
        )

        variable_sites = []

        alignment_length = alignment.get_alignment_length()

        for position in range(alignment_length):

            column = alignment[:, position]

            bases = set(column)

            bases.discard("-")
            bases.discard("N")
            bases.discard("?")

            if len(bases) > 1:

                variable_sites.append(position)

        return alignment, variable_sites

    def find_species_specific_snps(
        self,
        alignment_file: Path,
        reference_name: str,
    ):

        alignment, variable_sites = self.find_variable_sites(
            alignment_file
        )

        reference = None

        others = []

        for record in alignment:

            if reference_name.lower() in record.description.lower():

                reference = record

            else:

                others.append(record)

        if reference is None:

            raise ValueError(
                f"{reference_name} not found in alignment."
            )

        specific = []

        for position in variable_sites:

            ref_base = reference.seq[position]

            if ref_base in "-N?":
                continue

            unique = True

            for record in others:

                base = record.seq[position]

                if base == ref_base:

                    unique = False
                    break

            if unique:

                specific.append(
                    {
                        "position": position + 1,
                        "reference": ref_base,
                    }
                )

        return specific

    def summary(
        self,
        alignment_file: Path,
        reference_name: str,
    ):

        snps = self.find_species_specific_snps(
            alignment_file,
            reference_name,
        )

        print()

        print(f"Reference : {reference_name}")

        print(f"Species-specific SNPs : {len(snps)}")

        print()

        for snp in snps[:20]:

            print(snp)
