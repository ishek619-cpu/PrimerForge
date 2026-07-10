"""
BLAST output parser.
"""

from pathlib import Path

from primerforge.specificity.models import BlastHit


class BlastParser:
    """
    Parse BLAST tabular output into BlastHit objects.
    """

    def parse(
        self,
        blast_output: Path,
    ) -> list[BlastHit]:

        hits = []

        if not blast_output.exists():

            return hits

        with open(blast_output) as handle:

            for line in handle:

                if not line.strip():

                    continue

                fields = line.rstrip().split("\t")

                hits.append(

                    BlastHit(

                        accession=fields[1],

                        species=fields[1],

                        identity=float(fields[2]),

                        coverage=100.0,

                        alignment_length=int(fields[3]),

                        mismatches=int(fields[4]),

                        gap_opens=int(fields[5]),

                        qstart=int(fields[6]),

                        qend=int(fields[7]),

                        sstart=int(fields[8]),

                        send=int(fields[9]),

                        strand=fields[12],

                        bitscore=float(fields[11]),

                        evalue=float(fields[10]),

                    )

                )

        return hits
