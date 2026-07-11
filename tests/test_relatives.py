from primerforge.species.taxonomy import TaxonomyRecord
from primerforge.species.relatives import RelativeSpeciesFinder


def test_genus():

    finder = RelativeSpeciesFinder(
        email="test@test.com",
    )

    record = TaxonomyRecord(

        taxid="1",

        scientific_name="Oreochromis niloticus",

        rank="species",

        lineage=[
            "Eukaryota",
            "Chordata",
            "Oreochromis",
        ],

    )

    assert finder.genus(record) == "Oreochromis"
