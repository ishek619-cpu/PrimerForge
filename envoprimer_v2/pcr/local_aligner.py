"""
Sliding local aligner for PCR primer binding.

The aligner searches every possible placement of a primer within a
candidate window and returns the highest-scoring binding site.
"""

from __future__ import annotations

from envoprimer_v2.pcr.models import BindingSite


class LocalAligner:
    """
    Sliding local aligner.

    The primer is aligned against every possible position inside the
    candidate window. The best alignment is returned.
    """

    def __init__(self, terminal_region: int = 5):
        self.terminal_region = terminal_region

    def align(
        self,
        *,
        primer: str,
        window: str,
        sequence_id: int,
        strand: str,
        window_start: int,
    ) -> BindingSite:

        primer = primer.upper()
        window = window.upper()

        p_len = len(primer)

        if len(window) < p_len:
            raise ValueError("Window shorter than primer.")

        best = None
        best_score = -1.0

        for offset in range(len(window) - p_len + 1):

            template = window[offset:offset + p_len]

            mismatches = 0
            terminal = 0

            for i, (p, t) in enumerate(zip(primer, template)):

                if p != t:
                    mismatches += 1

                    if i >= p_len - self.terminal_region:
                        terminal += 1

            identity = 100.0 * (p_len - mismatches) / p_len

            score = (
                identity
                - (terminal * 10)
                - ((mismatches - terminal) * 2)
            )

            if score > best_score:

                best_score = score

                best = BindingSite(
                    sequence_id=sequence_id,
                    strand=strand,
                    start=window_start + offset,
                    end=window_start + offset + p_len,
                    identity=identity,
                    mismatches=mismatches,
                    terminal_mismatches=terminal,
                    score=score,
                    primer_alignment=primer,
                    template_alignment=template,
                )

        return best
