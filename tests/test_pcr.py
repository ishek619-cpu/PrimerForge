"""
Tests for PCR product prediction.
"""

from envoprimer.specificity.pcr import (
    PCRProductFinder,
)


def test_pcr_product_prediction():

    forward_hits = [

        {

            "subject": "chr1",

            "subject_start": 100,

            "subject_end": 120,

            "query_start": 1,

            "query_end": 20,

            "identity": 100.0,

            "strand": "plus",

        }

    ]

    reverse_hits = [

        {

            "subject": "chr1",

            "subject_start": 250,

            "subject_end": 270,

            "query_start": 1,

            "query_end": 20,

            "identity": 99.0,

            "strand": "minus",

        }

    ]

    finder = PCRProductFinder()

    products = finder.find_products(

        forward_hits,

        reverse_hits,

        min_size=50,

        max_size=500,

    )

    assert len(products) == 1

    product = products[0]

    assert product.chromosome == "chr1"

    assert product.forward_start == 100

    assert product.forward_end == 120

    assert product.reverse_start == 250

    assert product.reverse_end == 270

    assert product.size == 171

    assert product.identity == 99.0

    assert product.strand == "plus"
