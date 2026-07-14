"""
PCR primer alignment.

Performs ungapped primer-to-template alignment and reports
binding quality for downstream PCR evaluation.
"""

from __future__ import annotations

from envoprimer.pcr.models import PrimerBinding


class PrimerAligner:

    def reverse_complement(
        self,
        sequence: str,
    ) -> str:

        table = str.maketrans(
            "ACGTN",
            "TGCAN",
        )

        return sequence.translate(table)[::-1]

    ########################################################

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

        if strand == "-":
            primer = self.reverse_complement(
                primer,
            )

        end = start + len(primer)

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

        region = template[start:end]

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

        ####################################################
        # Count consecutive mismatches from the primer 3′ end
        ####################################################

        terminal = 0

        for p, t in zip(
            reversed(primer),
            reversed(region),
        ):

            if p == t:
                break

            terminal += 1

        identity = (
            matches / len(primer)
        ) * 100.0

        #
        # Simple binding score for now.
        # We'll replace this later with
        # nearest-neighbor thermodynamics.
        #

        score = identity
        score -= mismatches * 10
        score -= terminal * 20

        if score < 0:
            score = 0.0

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
