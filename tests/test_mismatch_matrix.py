"""
Tests for mismatch thermodynamic penalty matrix.
"""

from primerforge.specificity.mismatch_matrix import (
    MATCH_SCORE,
    mismatch_penalty,
)


def test_perfect_matches():

    assert mismatch_penalty("A", "A") == MATCH_SCORE

    assert mismatch_penalty("C", "C") == MATCH_SCORE

    assert mismatch_penalty("G", "G") == MATCH_SCORE

    assert mismatch_penalty("T", "T") == MATCH_SCORE


def test_wobble_pair():

    assert mismatch_penalty("G", "T") == 0.25

    assert mismatch_penalty("T", "G") == 0.25


def test_moderate_penalties():

    assert mismatch_penalty("A", "C") == 0.50

    assert mismatch_penalty("C", "A") == 0.50

    assert mismatch_penalty("A", "G") == 0.60

    assert mismatch_penalty("G", "A") == 0.60


def test_severe_penalties():

    assert mismatch_penalty("A", "T") == 0.90

    assert mismatch_penalty("T", "A") == 0.90

    assert mismatch_penalty("C", "G") == 0.90

    assert mismatch_penalty("G", "C") == 0.90


def test_unknown_bases():

    assert mismatch_penalty("N", "A") == 1.0

    assert mismatch_penalty("A", "N") == 1.0

    assert mismatch_penalty("N", "N") == 1.0
