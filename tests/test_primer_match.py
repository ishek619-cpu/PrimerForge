from envoprimer.analysis.primer_match import (
    PrimerMatcher,
)


def test_perfect_match():

    matcher = PrimerMatcher()

    result = matcher.match(
        "ATGC",
        "ATGCATGC",
        0,
    )

    assert result.matched
    assert result.identity == 100.0
    assert result.mismatches == 0


def test_gap_handling():

    matcher = PrimerMatcher()

    result = matcher.match(
        "ATGC",
        "AT-GCATGC",
        0,
    )

    assert result.matched
    assert result.identity == 100.0


def test_single_mismatch():

    matcher = PrimerMatcher()

    result = matcher.match(
        "ATGC",
        "ATGTATGC",
        0,
    )

    assert result.matched
    assert result.mismatches == 1


def test_multiple_mismatches():

    matcher = PrimerMatcher()

    result = matcher.match(
        "ATGC",
        "TTTTATGC",
        0,
    )

    assert not result.matched
    assert result.mismatches > 1


def test_short_sequence():

    matcher = PrimerMatcher()

    result = matcher.match(
        "ATGCATGC",
        "ATGC",
        0,
    )

    assert not result.matched
