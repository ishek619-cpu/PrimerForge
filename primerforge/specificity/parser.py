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

        with open(
            blast_output,
            encoding="utf-8",
        ) as handle:

            for line in handle:

                line = line.strip()

                if not line:
                    continue

                fields = line.split("\t")

                if len(fields) < 15:
                    continue

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

                        #
                        # NEW
                        #
                        "query_sequence": fields[10],

                        "subject_sequence": fields[11],

                        "evalue": float(fields[12]),

                        "bitscore": float(fields[13]),

                        "strand": fields[14],

                    }
                )

        return hits
