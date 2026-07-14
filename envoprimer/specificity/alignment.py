"""
Primer alignment scoring.

Provides both the new indexed alignment engine and
backward compatibility for the legacy ThreePrimeAnalyzer.
"""

from __future__ import annotations

from dataclasses import dataclass


# ======================================================================
# New alignment engine
# ======================================================================

@dataclass(slots=True)
class AlignmentResult:

    matches: int

    mismatches: int

    terminal_mismatches: int

    consecutive_terminal_mismatches: int

    identity: float

    score: float

    passed: bool


class PrimerAlignment:

    def __init__(

        self,

        max_mismatches: int = 2,

        min_identity: float = 90.0,

    ):

        self.max_mismatches = max_mismatches

        self.min_identity = min_identity

    ############################################################

    def reverse_complement(

        self,

        sequence: str,

    ) -> str:

        table = str.maketrans(

            "ACGTN",

            "TGCAN",

        )

        return sequence.translate(table)[::-1]

    ############################################################

    def align(

        self,

        primer: str,

        target: str,

        start: int,

        strand: str = "+",

    ) -> AlignmentResult:

        primer = primer.upper()

        target = target.upper()

        if strand == "-":

            primer = self.reverse_complement(

                primer,

            )

        end = start + len(primer)

        if start < 0 or end > len(target):

            return AlignmentResult(

                matches=0,

                mismatches=len(primer),

                terminal_mismatches=len(primer),

                consecutive_terminal_mismatches=len(primer),

                identity=0.0,

                score=0.0,

                passed=False,

            )

        region = target[start:end]

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

        terminal = 0

        consecutive = 0

        for p, t in zip(

            reversed(primer),

            reversed(region),

        ):

            if p != t:

                terminal += 1

                consecutive += 1

            else:

                break

        ########################################################

        identity = (

            matches / len(primer)

        ) * 100.0

        score = identity

        score -= mismatches * 8

        score -= terminal * 15

        score -= consecutive * 10

        score = max(

            0.0,

            score,

        )

        passed = (

            mismatches <= self.max_mismatches

            and identity >= self.min_identity

        )

        return AlignmentResult(

            matches=matches,

            mismatches=mismatches,

            terminal_mismatches=terminal,

            consecutive_terminal_mismatches=consecutive,

            identity=round(identity, 2),

            score=round(score, 2),

            passed=passed,

        )


# ======================================================================
# Legacy compatibility
# ======================================================================

@dataclass(slots=True)
class AlignmentPosition:

    primer_position: int

    primer_base: str

    subject_base: str

    aligned: bool

    mismatch: bool


class AlignmentReconstructor:
    """
    Compatibility layer for the old ThreePrimeAnalyzer.
    """

    def reconstruct(

        self,

        hit,

        primer_length: int,

    ):

        query = (

            hit.query_sequence.upper()

            if hit.query_sequence

            else ""

        )

        subject = (

            hit.subject_sequence.upper()

            if hit.subject_sequence

            else ""

        )

        alignment = []

        for i in range(primer_length):

            if i < len(query):

                q = query[i]

            else:

                q = "-"

            if i < len(subject):

                s = subject[i]

            else:

                s = "-"

            aligned = (

                q != "-"

                and s != "-"

            )

            alignment.append(

                AlignmentPosition(

                    primer_position=i,

                    primer_base=q,

                    subject_base=s,

                    aligned=aligned,

                    mismatch=(

                        aligned

                        and q != s

                    ),

                )

            )

        return alignment
