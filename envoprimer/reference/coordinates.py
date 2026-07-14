"""
Universal coordinate system utilities.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class Coordinate:

    start: int

    end: int

    system: str = "reference"

    @property
    def length(self) -> int:

        return self.end - self.start + 1

    def contains(
        self,
        position: int,
    ) -> bool:

        return self.start <= position <= self.end

    def overlaps(
        self,
        other: "Coordinate",
    ) -> bool:

        if self.system != other.system:

            return False

        return not (
            self.end < other.start
            or
            other.end < self.start
        )

    def shift(
        self,
        offset: int,
    ) -> "Coordinate":

        return Coordinate(
            start=self.start + offset,
            end=self.end + offset,
            system=self.system,
        )

    def to_dict(self) -> dict:

        return {

            "start": self.start,

            "end": self.end,

            "length": self.length,

            "system": self.system,

        }

    def __repr__(self):

        return (
            f"Coordinate("
            f"{self.start}-{self.end}, "
            f"system='{self.system}')"
        )
