"""
Tests for alignment reconstruction.
"""

from envoprimer.specificity.alignment import AlignmentReconstructor
from envoprimer.specificity.models import BlastHit


def test_alignment_reconstruction():

    hit = BlastHit(

        accession="ABC123",

        species="Oreochromis niloticus",

        identity=95.0,

        coverage=100.0,

        alignment_length=18,

        mismatches=1,

        gap_opens=0,

        qstart=2,

        qend=19,

        sstart=100,

        send=117,

        query_sequence="TCGATCGATCGATCGATC",

        subject_sequence="TCGATCGATCGATGGATC",

        strand="plus",

        bitscore=42.0,

        evalue=1e-20,

    )

    reconstructor = AlignmentReconstructor()

    positions = reconstructor.reconstruct(
        hit,
        primer_length=20,
    )

    assert len(positions) == 20

    #
    # First base is unaligned
    #
    assert positions[0].aligned is False

    #
    # Last base is unaligned
    #
    assert positions[19].aligned is False

    #
    # Alignment starts at primer position 2
    #
    assert positions[1].aligned is True

    #
    # Exactly one mismatch
    #
    mismatches = [
        p for p in positions
        if p.mismatch
    ]

    assert len(mismatches) == 1
