"""
Primer pair deduplication.
"""

from envoprimer.models.pair import PrimerPair


class PrimerPairDeduplicator:
    """
    Remove duplicate primer pairs.
    """

    def deduplicate(
        self,
        pairs: list[PrimerPair],
    ) -> list[PrimerPair]:

        unique = {}

        for pair in pairs:

            key = (
                pair.forward.sequence,
                pair.reverse.sequence,
            )

            if key not in unique:
                unique[key] = pair

        return list(unique.values())
