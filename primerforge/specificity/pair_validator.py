"""
Pair-level PCR specificity validation.
"""

from primerforge.models.pair import PrimerPair
from primerforge.specificity.models import SpecificityResult


class PairSpecificityValidator:
    """
    Validate whether a primer pair can amplify
    the same off-target sequence.
    """

    def validate(
        self,
        pair: PrimerPair,
    ) -> dict:

        result = pair.specificity_result

        #
        # No specificity information available.
        #
        if result is None:

            return {

                "passed": True,

                "off_target_amplicons": 0,

                "closest_off_target": None,

                "predicted_product_size": None,

            }

        #
        # Collect off-target species.
        #

        forward_species = {

            hit.species

            for hit in result.forward_hits

        }

        reverse_species = {

            hit.species

            for hit in result.reverse_hits

        }

        #
        # Species hit by BOTH primers.
        #

        shared = forward_species.intersection(
            reverse_species,
        )

        if not shared:

            return {

                "passed": True,

                "off_target_amplicons": 0,

                "closest_off_target": None,

                "predicted_product_size": None,

            }

        #
        # Potential off-target amplification.
        #

        species = sorted(
            shared,
        )[0]

        return {

            "passed": False,

            "off_target_amplicons": len(shared),

            "closest_off_target": species,

            "predicted_product_size": pair.product_size,

        }
