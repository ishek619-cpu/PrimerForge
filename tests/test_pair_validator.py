"""
Tests for pair-level specificity validator.
"""

from primerforge.models.pair import PrimerPair
from primerforge.models.primer import Primer

from primerforge.specificity.models import (
    BlastHit,
    SpecificityResult,
)

from primerforge.specificity.pair_validator import (
    PairSpecificityValidator,
)


def test_pair_validator():

    pair = PrimerPair(

        forward=Primer(
            sequence="AAAAAAAAAAAAAAAAAAAA",
            start=1,
            end=20,
            strand="+",
            length=20,
        ),

        reverse=Primer(
            sequence="TTTTTTTTTTTTTTTTTTTT",
            start=120,
            end=140,
            strand="-",
            length=20,
        ),

        product_size=120,
    )

    hit = BlastHit(

        accession="NC_000001",

        species="Oreochromis aureus",

        identity=100.0,

        coverage=100.0,

        alignment_length=20,

        mismatches=0,

        gap_opens=0,

        qstart=1,

        qend=20,

        sstart=100,

        send=119,

        strand="+",

        bitscore=40,

        evalue=0.0,
    )

    pair.specificity_result = SpecificityResult(

        forward_hits=[hit],

        reverse_hits=[hit],
    )

    validator = PairSpecificityValidator()

    result = validator.validate(
        pair,
    )

    assert result["passed"] is False

    assert result["off_target_amplicons"] == 1

    assert result["closest_off_target"] == "Oreochromis aureus"
