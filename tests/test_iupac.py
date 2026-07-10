from primerforge.analysis.iupac import (
    expand,
    matches,
    mismatch,
    identity,
    mismatch_count,
)


def test_expand():

    assert expand("A") == {"A"}

    assert expand("R") == {"A", "G"}

    assert expand("N") == {
        "A",
        "C",
        "G",
        "T",
    }


def test_matches():

    assert matches("A", "A")

    assert matches("R", "A")

    assert matches("R", "G")

    assert matches("N", "T")

    assert not matches("A", "C")


def test_identity():

    assert identity(
        "ATGC",
        "ATGC",
    ) == 100.0

    assert identity(
        "ATGC",
        "ATGT",
    ) == 75.0

    assert identity(
        "ATGR",
        "ATGA",
    ) == 100.0


def test_mismatch_count():

    assert mismatch_count(
        "ATGC",
        "ATGC",
    ) == 0

    assert mismatch_count(
        "ATGC",
        "ATGT",
    ) == 1

    assert mismatch_count(
        "ATGR",
        "ATGA",
    ) == 0


def test_unknown_base():

    assert matches(
        "X",
        "A",
    )

    assert identity(
        "AXGC",
        "ATGC",
    ) == 100.0
