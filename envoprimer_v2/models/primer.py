"""
Primer model.

Used throughout EnvoPrimer.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class Primer:

    sequence: str

    start: int

    end: int

    strand: str

    tm: float

    gc: float

    length: int

    penalty: float = 0.0

    hairpin_tm: float = 0.0

    homodimer_tm: float = 0.0

    heterodimer_tm: float = 0.0

    population_score: float = 0.0

    specificity_score: float = 0.0

    three_prime_score: float = 0.0

    thermo_score: float = 0.0

    final_score: float = 0.0

    notes: list[str] = field(default_factory=list)

    @property
    def three_prime_position(self) -> int:

        if self.strand == "+":

            return self.end

        return self.start

    def to_dict(self):

        return {

            "sequence": self.sequence,

            "start": self.start,

            "end": self.end,

            "strand": self.strand,

            "tm": self.tm,

            "gc": self.gc,

            "length": self.length,

            "penalty": self.penalty,

            "population_score": self.population_score,

            "specificity_score": self.specificity_score,

            "three_prime_score": self.three_prime_score,

            "thermo_score": self.thermo_score,

            "final_score": self.final_score,

        }
