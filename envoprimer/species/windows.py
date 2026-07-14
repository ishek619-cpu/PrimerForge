"""
Diagnostic window generation.

Builds candidate sequence windows around diagnostic SNPs for
species-specific primer design.
"""

from __future__ import annotations

from dataclasses import dataclass

from envoprimer.models.region import Region


@dataclass(slots=True)
class DiagnosticWindow:
    """
    Candidate window for Primer3.
    """

    start: int
    end: int
    sequence: str
    diagnostic_sites: list[int]

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_region(self) -> Region:
        """
        Convert a diagnostic window into a Primer3 Region.
        """

        return Region(
            start=self.start,
            end=self.end,
            length=self.length,
            score=float(len(self.diagnostic_sites)),
        )


class DiagnosticWindowBuilder:
    """
    Build sequence windows centred on diagnostic SNPs.
    """

    def __init__(
        self,
        window_size: int = 180,
        overlap: int = 40,
    ):

        self.window_size = window_size
        self.overlap = overlap

    def build(
        self,
        reference_sequence: str,
        diagnostic_sites,
    ) -> list[DiagnosticWindow]:

        if not diagnostic_sites:
            return []

        length = len(reference_sequence)

        positions = sorted(
            site.position
            for site in diagnostic_sites
        )

        windows: list[DiagnosticWindow] = []

        used = set()

        for centre in positions:

            start = max(
                0,
                centre - self.window_size // 2,
            )

            end = min(
                length,
                start + self.window_size,
            )

            start = max(
                0,
                end - self.window_size,
            )

            key = (
                start,
                end,
            )

            if key in used:
                continue

            used.add(key)

            contained = [

                pos

                for pos in positions

                if start <= pos < end

            ]

            windows.append(

                DiagnosticWindow(

                    start=start,

                    end=end,

                    sequence=reference_sequence[start:end],

                    diagnostic_sites=contained,

                )

            )

        windows.sort(

            key=lambda w: (

                -len(w.diagnostic_sites),

                w.start,

            )

        )

        merged: list[DiagnosticWindow] = []

        for window in windows:

            if not merged:

                merged.append(window)

                continue

            previous = merged[-1]

            overlap = previous.end - window.start

            if overlap >= self.overlap:

                continue

            merged.append(window)

        return merged
