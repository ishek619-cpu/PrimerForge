"""
Base rule definitions.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class RuleResult:

    name: str

    passed: bool

    score: float

    reason: str = ""


class Rule:

    name = "Unnamed Rule"

    mandatory = True

    def evaluate(self, pair):

        raise NotImplementedError
