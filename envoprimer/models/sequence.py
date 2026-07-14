"""
EnvoPrimer Sequence Model
"""

from dataclasses import dataclass


@dataclass(slots=True)
class SequenceRecord:
    """
    Metadata describing a downloaded sequence.
    """

    accession: str
    organism: str
    gene: str
    length: int
    topology: str
    molecule: str
