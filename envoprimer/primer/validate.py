"""
Primer validation and filtering.
"""

from envoprimer.models.primer import Primer


class PrimerValidator:

    def __init__(self):
        pass

    def passes(self, primer: Primer) -> bool:

        if primer.length < 20:
            return False

        if primer.length > 24:
            return False

        if primer.tm < 58:
            return False

        if primer.tm > 62:
            return False

        if primer.gc < 40:
            return False

        if primer.gc > 60:
            return False

        if not primer.gc_clamp:
            return False

        if primer.homopolymer:
            return False

        return True

    def filter(
        self,
        primers: list[Primer],
    ) -> list[Primer]:

        accepted = []

        for primer in primers:

            if self.passes(primer):

                accepted.append(primer)

        return accepted
