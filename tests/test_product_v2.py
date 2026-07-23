from envoprimer_v2.pcr.models import BindingSite
from envoprimer_v2.pcr.product import ProductPredictor


forward = BindingSite(
    sequence_id=0,
    strand="+",
    start=100,
    end=120,
    identity=100,
    mismatches=0,
    terminal_mismatches=0,
    score=100,
    primer_alignment="",
    template_alignment="",
)

reverse = BindingSite(
    sequence_id=0,
    strand="-",
    start=250,
    end=270,
    identity=100,
    mismatches=0,
    terminal_mismatches=0,
    score=100,
    primer_alignment="",
    template_alignment="",
)

predictor = ProductPredictor()

products = predictor.predict([forward], [reverse])

print()
print("=" * 60)
print("PCR PRODUCT TEST")
print("=" * 60)
print("Products :", len(products))

if products:
    p = products[0]
    print("Amplicon :", p.start, "-", p.end)
    print("Size     :", p.size)
