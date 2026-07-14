"""
Ecological validation models.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class EcologicalHit:
    species: str
    taxonomy: str
    region: str
    source: str
    amplified: bool
    risk: float


@dataclass(slots=True)
class EcologicalResult:
    target_species: str
    total_species: int
    amplified_species: int
    ecological_score: float
    hits: list[EcologicalHit] = field(default_factory=list)
