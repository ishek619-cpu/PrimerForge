"""
BLAST output parser.
"""

from __future__ import annotations

from pathlib import Path

from envoprimer.specificity.models import BlastHit


class BlastParser:
    """
    Parse BLAST tabular output.
    """

    def parse(
        self,
        blast_output: Path,
    ) -> list[BlastHit]:

        hits: list[BlastHit] = []

        if not blast_output.exists():

            return hits

        with blast_output.open(
            encoding="utf-8",
        ) as handle:

            for line in handle:

                line = line.strip()

                if not line:

                    continue

                fields = line.split("\t")

                if len(fields) < 15:

                    continue

                sstart = int(fields[8])

                send = int(fields[9])

                strand = fields[14]

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

                        sstart=sstart,

                        send=send,

                        strand=strand,

                        bitscore=float(fields[13]),

                        evalue=float(fields[12]),

                        three_prime_mismatches=0,

                    )

                )

        return hits
