"""
Binding evaluation rules for PCR primer binding.
"""

from __future__ import annotations

from envoprimer_v2.pcr.models import BindingSite


class BindingEvaluator:
    """
    Decide whether a primer binding site is acceptable for PCR.
    """

    def __init__(
        self,
        min_identity: float = 90.0,
        max_mismatches: int = 2,
        max_terminal_mismatches: int = 1,
    ):
        self.min_identity = min_identity
        self.max_mismatches = max_mismatches
        self.max_terminal_mismatches = max_terminal_mismatches

    def accept(self, site: BindingSite) -> bool:

        if site.identity < self.min_identity:
            return False

        if site.mismatches > self.max_mismatches:
            return False

        if site.terminal_mismatches > self.max_terminal_mismatches:
            return False

        return True
