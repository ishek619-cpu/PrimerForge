from dataclasses import FrozenInstanceError

import pytest

from primerforge.specificity.models import (
    BlastHit,
    OffTargetHit,
    SpecificityResult,
)


def test_blast_hit_creation():

    hit = BlastHit(
        accession="NC_000001",
        species="Oreochromis niloticus",
        identity=100.0,
        coverage=100.0,
        alignment_length=20,
        mismatches=0,
        gap_opens=0,
        qstart=1,
        qend=20,
        sstart=100,
        send=119,
        strand="plus",
        bitscore=40.0,
        evalue=1e-10,
    )

    assert hit.species == "Oreochromis niloticus"
    assert hit.identity == 100.0


def test_offtarget_hit():

    hit = OffTargetHit(
        accession="NC_000002",
        species="Oreochromis mossambicus",
        identity=97.5,
        coverage=100.0,
        alignment_length=20,
        mismatches=1,
        gap_opens=0,
        qstart=1,
        qend=20,
        sstart=50,
        send=69,
        strand="plus",
        bitscore=36.0,
        evalue=1e-6,
        penalty=5.0,
    )

    assert hit.penalty == 5.0


def test_specificity_result_defaults():

    result = SpecificityResult()

    assert result.forward_hits == []
    assert result.reverse_hits == []
    assert result.target_hits == []
    assert result.off_target_hits == []
    assert result.specificity_score == 100.0
    assert result.passed is True


def test_models_are_frozen():

    hit = BlastHit(
        accession="A",
        species="Species",
        identity=100,
        coverage=100,
        alignment_length=20,
        mismatches=0,
        gap_opens=0,
        qstart=1,
        qend=20,
        sstart=1,
        send=20,
        strand="plus",
        bitscore=40,
        evalue=0,
    )

    with pytest.raises(FrozenInstanceError):
        hit.identity = 90
