"""
BLAST output parser.
"""

from pathlib import Path


class BlastParser:
    """
    Parse BLAST tabular output.
    """

    def parse(
        self,
        blast_output: Path,
    ):

        hits = []

        if not blast_output.exists():
            return hits

        with open(blast_output) as handle:

            for line in handle:

                if not line.strip():
                    continue

                fields = line.rstrip().split("\t")

                hits.append(
                    {

                        "query": fields[0],

                        "subject": fields[1],

                        "identity": float(fields[2]),

                        "length": int(fields[3]),

                        "mismatches": int(fields[4]),

                        "gapopen": int(fields[5]),

                        "query_start": int(fields[6]),

                        "query_end": int(fields[7]),

                        "subject_start": int(fields[8]),

                        "subject_end": int(fields[9]),

                        "evalue": float(fields[10]),

                        "bitscore": float(fields[11]),

                        "strand": fields[12],

                    }
                )

        return hits
