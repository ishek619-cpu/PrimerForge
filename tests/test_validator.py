from pathlib import Path

from primerforge.models.primer import Primer
from primerforge.models.pair import PrimerPair
from primerforge.analysis.snps import SNPFinder
from primerforge.validation.validator import PrimerValidator


def test_validator():

    pair = PrimerPair(

        forward=Primer(
            sequence="CCCTTCATCATTGCAGCTGC",
            start=556,
            end=575,
            strand="+",
            length=20,
        ),

        reverse=Primer(
            sequence="CTGAGTTTAGGCCTGTGGGG",
            start=640,
            end=659,
            strand="-",
            length=20,
        ),

        product_size=104,
    )

    snps = SNPFinder().find(
        Path(
            "data/alignments/alignment.fasta"
        )
    )

    result = PrimerValidator().validate(
        pair,
        Path(
            "data/alignments/alignment.fasta"
        ),
        snps,
    )

    assert result["final_score"] > 0
