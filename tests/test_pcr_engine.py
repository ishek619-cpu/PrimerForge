"""
Basic test for the PCR Engine.
"""

from envoprimer.pcr.engine import PCREngine


class Primer:
    def __init__(self, sequence):
        self.sequence = sequence


class PrimerPair:
    def __init__(self, forward, reverse):
        self.forward = Primer(forward)
        self.reverse = Primer(reverse)


#
# Background sequence
#
background = [
    "AAATGGAGGCTTTGGAAACTGAC"
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    "GACATGGCCTTCCCTCGAAT"
]

engine = PCREngine(
    minimum_product=20,
    maximum_product=200,
)

engine.build(background)

pair = PrimerPair(
    "TGGAGGCTTTGGAAACTGAC",
    "ATTCGAGGGAAGGCCATGTC",
)

products = engine.predict(pair)

print()

print("=" * 60)
print("PCR ENGINE TEST")
print("=" * 60)

print(f"Products detected : {len(products)}")

for p in products:
    print(
        f"Sequence {p.sequence_id} "
        f"Product {p.length} bp"
    )

