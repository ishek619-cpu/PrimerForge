from envoprimer_v2.pcr.engine import PCREngine
from envoprimer_v2.pcr.search import reverse_complement

# ---------------------------------------------------------------------
# Primers (always written 5'->3')
# ---------------------------------------------------------------------

forward = "TGGAGGCTTTGGAAACTGAC"

reverse = reverse_complement(forward)

# ---------------------------------------------------------------------
# Genome
#
# IMPORTANT:
# The genome contains the forward primer sequence and the
# REVERSE-COMPLEMENT of the reverse primer (i.e. the template strand).
# ---------------------------------------------------------------------

sequence = (
    "AAA"
    + forward
    + ("A" * 120)
    + reverse_complement(reverse)
    + "AAA"
)

print("=" * 70)
print("PCR ENGINE V2")
print("=" * 70)

print("Forward primer :", forward)
print("Reverse primer :", reverse)

print("\nGenome tail:")
print(sequence[-40:])

engine = PCREngine()
engine.build([sequence])

print("\nSearching forward sites...")
forward_sites = engine._find_sites(forward, "+")
print("Forward sites :", len(forward_sites))

for s in forward_sites:
    print(s)

print("\nSearching reverse sites...")
reverse_sites = engine._find_sites(reverse, "-")
print("Reverse sites :", len(reverse_sites))

for s in reverse_sites:
    print(s)

print("\nPredicting PCR products...")
products = engine.predict(forward, reverse)

print("Products :", len(products))

for i, p in enumerate(products, 1):
    print()
    print(f"Product {i}")
    print("Start :", p.start)
    print("End   :", p.end)
    print("Size  :", p.size)
