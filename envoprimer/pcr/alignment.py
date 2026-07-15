"""
PCR primer alignment.

Performs strand-aware ungapped alignment of primers to DNA templates.
"""

from __future__ import annotations

from envoprimer.pcr.models import PrimerBinding


class PrimerAligner:

    ############################################################

    def reverse_complement(
        self,
        sequence: str,
    ) -> str:

        table = str.maketrans(
            "ACGTN",
            "TGCAN",
        )

        return sequence.upper().translate(table)[::-1]

    ############################################################

    def align(
        self,
        primer: str,
        template: str,
        sequence_id: int,
        start: int,
        strand: str,
    ) -> PrimerBinding:

        primer = primer.upper()
        template = template.upper()

        end = start + len(primer)

        ########################################################

        if start < 0 or end > len(template):

            return PrimerBinding(
                sequence_id=sequence_id,
                strand=strand,
                start=start,
                end=end,
                mismatches=len(primer),
                terminal_mismatches=len(primer),
                identity=0.0,
                score=0.0,
            )

        ########################################################
        # Extract template region
        ########################################################

        region = template[start:end]

        ########################################################
        # IMPORTANT
        #
        # Reverse primer binds to the reverse-complement
        # of the template region.
        ########################################################

        if strand == "-":

            region = self.reverse_complement(
                region,
            )

        ########################################################

        matches = 0
        mismatches = 0

        for p, t in zip(

            primer,

            region,

        ):

            if p == t:

                matches += 1

            else:

                mismatches += 1

        ########################################################
        # Count consecutive mismatches from primer 3'
        ########################################################

        terminal = 0

        for p, t in zip(

            reversed(primer),

            reversed(region),

        ):

            if p == t:

                break

            terminal += 1

        ########################################################

        identity = (

            matches

            / len(primer)

        ) * 100.0

        ########################################################

        score = identity

        score -= mismatches * 10

        score -= terminal * 20

        if score < 0:

            score = 0.0

        ########################################################

        return PrimerBinding(

            sequence_id=sequence_id,

            strand=strand,

            start=start,

            end=end,

            mismatches=mismatches,

            terminal_mismatches=terminal,

            identity=round(identity, 2),

            score=round(score, 2),

        )
