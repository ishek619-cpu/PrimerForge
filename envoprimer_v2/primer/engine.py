"""
EnvoPrimer Engine.

Coordinates Primer3 and primer pair generation.

Population conservation, specificity and final scoring are now
performed by the new EnvoPrimer engines.
"""

from __future__ import annotations

from collections import Counter

from envoprimer_v2.primer3.designer import PrimerDesigner
from envoprimer_v2.primer3.parser import PrimerParser
from envoprimer_v2.primer.generator import PrimerPairGenerator

from envoprimer_v2.rules.runner import RuleEngine

from envoprimer_v2.rules.product_size import ProductSizeRule
from envoprimer_v2.rules.tm import TmRule
from envoprimer_v2.rules.gc import GCRule
from envoprimer_v2.rules.homopolymer import HomopolymerRule
from envoprimer_v2.rules.gc_clamp import GCClampRule
from envoprimer_v2.rules.poly3 import ThreePrimeRule


class EnvoPrimerEngine:

    def __init__(self):

        self.designer = PrimerDesigner()

        self.parser = PrimerParser()

        self.generator = PrimerPairGenerator()

        self.rules = RuleEngine()

        #
        # Only intrinsic primer quality rules remain here.
        #

        self.rules.add(ProductSizeRule())

        self.rules.add(TmRule())

        self.rules.add(GCRule())

        self.rules.add(HomopolymerRule())

        self.rules.add(GCClampRule())

        self.rules.add(ThreePrimeRule())

    def design(

        self,

        window,

        diagnostic_snps,

        background_sequences,

    ):

        print()
        print("=" * 60)
        print("Primer3 Design")
        print("=" * 60)

        result = self.designer.design(

            window,

        )

        forward, reverse = self.parser.parse(

            result,

            window.start,

        )

        print(f"Forward primers returned : {len(forward)}")
        print(f"Reverse primers returned : {len(reverse)}")

        pairs = self.generator.generate(

            forward,

            reverse,

        )

        print(f"Primer pairs generated   : {len(pairs)}")

        accepted = []

        rejected = Counter()

        for pair in pairs:

            outcome = self.rules.evaluate(

                pair,

            )

            if not outcome.passed:

                rejected[outcome.failed_rule] += 1

                continue

            #
            # Remaining scoring is handled later by:
            #
            #   ConservationEngine
            #   SpecificityEngine
            #   ScoringEngine
            #

            accepted.append(pair)

        print()

        print("=" * 60)
        print("Primer Evaluation Summary")
        print("=" * 60)

        print(f"Forward primers : {len(forward)}")
        print(f"Reverse primers : {len(reverse)}")
        print(f"Primer pairs    : {len(pairs)}")
        print(f"Accepted pairs  : {len(accepted)}")

        if rejected:

            print()

            print("Rejected by rule")

            for rule, count in rejected.items():

                print(f"  {rule:<25} {count}")

        print()

        return accepted
