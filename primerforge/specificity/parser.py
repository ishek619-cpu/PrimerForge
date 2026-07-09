"""
BLAST output parser.
"""

from pathlib import Path


class BlastParser:
    """
    Parse BLAST tabular output (outfmt 6).
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
                        "evalue": float(fields[10]),
                        "bitscore": float(fields[11]),
                    }
                )

        return hits
