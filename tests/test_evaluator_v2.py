from envoprimer_v2.pcr.models import BindingSite
from envoprimer_v2.pcr.evaluator import BindingEvaluator


site = BindingSite(
    sequence_id=0,
    strand="+",
    start=3,
    end=23,
    identity=100.0,
    mismatches=0,
    terminal_mismatches=0,
    score=100.0,
    primer_alignment="AAAA",
    template_alignment="AAAA",
)

evaluator = BindingEvaluator()

print()
print("=" * 60)
print("EVALUATOR TEST")
print("=" * 60)

print("Accepted:", evaluator.accept(site))
