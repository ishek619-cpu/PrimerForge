"""
Primer specificity analysis.
"""

from pathlib import Path

from envoprimer.models.primer import Primer

from envoprimer.specificity.blast import BlastRunner
from envoprimer.specificity.parser import BlastParser
from envoprimer.specificity.scorer import SpecificityScorer
from envoprimer.specificity.models import (
    OffTargetHit,
    SpecificityResult,
)


class SpecificityAnalyzer:
    """
    Analyse primer specificity using BLAST.
    """

    def __init__(self):

        self.blast = BlastRunner()

        self.parser = BlastParser()

        self.scorer = SpecificityScorer()

    def analyse(
        self,
        primer: Primer,
        database: str,
        blast_output: Path,
        target_species: str,
    ) -> SpecificityResult:

        blast_output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        #
        # Write primer FASTA for BLAST
        #
        query = blast_output.with_suffix(".fa")

        with query.open(
            "w",
            encoding="utf-8",
        ) as handle:

            handle.write(
                f">{primer.sequence}\n"
            )

            handle.write(
                f"{primer.sequence}\n"
            )

        #
        # Run BLAST
        #
        self.blast.search(

            query=query,

            database=database,

            output=blast_output,

        )

        hits = self.parser.parse(
            blast_output,
        )

        target_hits = []

        off_target_hits = []

        for hit in hits:

            if target_species.lower() in hit.species.lower():

                target_hits.append(
                    hit,
                )

            else:

                off_target_hits.append(

                    OffTargetHit(

                        accession=hit.accession,

                        species=hit.species,

                        identity=hit.identity,

                        coverage=hit.coverage,

                        alignment_length=hit.alignment_length,

                        mismatches=hit.mismatches,

                        gap_opens=hit.gap_opens,

                        qstart=hit.qstart,

                        qend=hit.qend,

                        sstart=hit.sstart,

                        send=hit.send,

                        strand=hit.strand,

                        bitscore=hit.bitscore,

                        evalue=hit.evalue,

                        three_prime_mismatches=hit.three_prime_mismatches,

                        penalty=100.0 - hit.identity,

                    )

                )

        specificity = self.scorer.score(
            hits,
        )

        return SpecificityResult(

            forward_hits=hits,

            reverse_hits=[],

            target_hits=target_hits,

            off_target_hits=off_target_hits,

            specificity_score=specificity,

            passed=specificity >= 90.0,

            rejection_reason=(
                None
                if specificity >= 90.0
                else "Low species specificity"
            ),

        )
