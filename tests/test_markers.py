from primerforge.species.markers import (
    marker_synonyms,
    build_query,
)


def test_marker_synonyms():

    synonyms = marker_synonyms(
        "COI",
    )

    assert "COI" in synonyms

    assert "COX1" in synonyms


def test_query_builder():

    query = build_query(
        "Oreochromis niloticus",
        "COI",
    )

    assert "Oreochromis niloticus" in query

    assert "COX1" in query

    assert "COI" in query
