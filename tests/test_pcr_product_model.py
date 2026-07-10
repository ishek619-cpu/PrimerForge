"""
Tests for PCRProduct model.
"""

from primerforge.models.pcr_product import (
    PCRProduct,
)


def test_pcr_product():

    product = PCRProduct(

        chromosome="NC_013663",

        forward_start=100,

        reverse_end=250,

        product_size=151,

        forward_strand="+",

        reverse_strand="-",

    )

    assert product.chromosome == "NC_013663"

    assert product.product_size == 151

    assert product.passed


def test_failed_product():

    product = PCRProduct(

        chromosome="chr1",

        forward_start=10,

        reverse_end=600,

        product_size=591,

        forward_strand="+",

        reverse_strand="-",

        passed=False,

        reason="Too large",

    )

    assert not product.passed

    assert product.reason == "Too large"
