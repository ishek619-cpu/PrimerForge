"""
PCR binding evaluator.

Determines whether a primer binding is capable
of supporting PCR amplification.
"""

from __future__ import annotations

from envoprimer.pcr.models import PrimerBinding


class PCREvaluator:

    def __init__(
        self,
        max_mismatches: int = 2,
        max_terminal_mismatches: int = 1,
        minimum_identity: float = 90.0,
    ):

        self.max_mismatches = max_mismatches
        self.max_terminal_mismatches = max_terminal_mismatches
        self.minimum_identity = minimum_identity

    ########################################################

    def passes(
        self,
        binding: PrimerBinding,
    ) -> bool:

        #
        # Too many mismatches
        #

        if binding.mismatches > self.max_mismatches:
            return False

        #
        # Too many 3' mismatches
        #

        if (
            binding.terminal_mismatches
            > self.max_terminal_mismatches
        ):
            return False

        #
        # Identity too low
        #

        if binding.identity < self.minimum_identity:
            return False

        return True

