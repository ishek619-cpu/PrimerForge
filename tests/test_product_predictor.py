"""
Tests for PCR product prediction.
"""

from primerforge.pcr.product_predictor import (
    PCRProductPredictor,
)
from primerforge.specificity.models import (
    BlastHit,
)


def make_hit(start, end):

    return BlastHit(

        accession="NC_013663",

        species="Oreochromis niloticus",

        identity=100.0,

        coverage=100.0,

        alignment_length=20,

        mismatches=0,

        gap_opens=0,

        qstart=1,

        qend=20,

        sstart=start,

        send=end,

        query_sequence="ATGCATGCATGCATGCATGC",

        subject_sequence="ATGCATGCATGCATGCATGC",

        strand="plus",

        bitscore=40.0,

        evalue=0.0,

    )


def test_predict():

    predictor = PCRProductPredictor()

    product = predictor.predict(

        make_hit(100, 119),

        make_hit(250, 269),

    )

    assert product is not None

    assert product.passed

    assert product.chromosome == "NC_013663"

    assert product.product_size == 170


def test_wrong_order():

    predictor = PCRProductPredictor()

    product = predictor.predict(

        make_hit(300, 319),

        make_hit(150, 169),

    )

    assert product is None


def test_large_product():

    predictor = PCRProductPredictor(

        max_product=200,

    )

    product = predictor.predict(

        make_hit(100, 119),

        make_hit(500, 519),

    )

    assert product is not None

    assert not product.passed

    assert product.reason == "Product too large"
