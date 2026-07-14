"""
PCR product prediction.
"""

from __future__ import annotations

from envoprimer.pcr.models import Amplicon, PrimerBinding


class PCRProductPredictor:

    def __init__(

        self,

        minimum_product: int = 70,

        maximum_product: int = 200,

    ):

        self.minimum_product = minimum_product

        self.maximum_product = maximum_product

    ############################################################

    def predict(

        self,

        forward: list[PrimerBinding],

        reverse: list[PrimerBinding],

    ) -> list[Amplicon]:

        products = []

        ########################################################

        reverse_lookup = {}

        for r in reverse:

            reverse_lookup.setdefault(

                r.sequence_id,

                [],

            ).append(r)

        ########################################################

        for f in forward:

            if f.sequence_id not in reverse_lookup:

                continue

            ####################################################
            # Forward primer must bind to the + strand
            ####################################################

            if f.strand != "+":

                continue

            for r in reverse_lookup[

                f.sequence_id

            ]:

                ################################################
                # Reverse primer must bind to the - strand
                ################################################

                if r.strand != "-":

                    continue

                ################################################
                # Reverse primer must be downstream
                ################################################

                if r.start <= f.start:

                    continue

                ################################################
                # Product length
                ################################################

                length = (

                    r.end

                    - f.start

                )

                if length < self.minimum_product:

                    continue

                if length > self.maximum_product:

                    continue

                ################################################

                products.append(

                    Amplicon(

                        sequence_id=f.sequence_id,

                        forward=f,

                        reverse=r,

                        length=length,

                        passed=True,

                    )

                )

        return products
