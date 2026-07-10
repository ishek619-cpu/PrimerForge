"""
Tests for PCR product prediction.
"""

from primerforge.specificity.pcr import (
    PCRProductFinder,
)


def test_pcr_product_prediction():

    forward_hits = [

        {

            "subject": "chr1",

            "query_start": 100,

            "query_end": 120,

            "identity": 100.0,

        }

    ]

    reverse_hits = [

        {

            "subject": "chr1",

            "query_start": 250,

            "query_end": 270,

            "identity": 99.0,

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

    assert products[0].size == 171

    assert products[0].chromosome == "chr1"

    assert products[0].identity == 99.0
