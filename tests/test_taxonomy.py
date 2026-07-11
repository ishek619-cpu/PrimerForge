from primerforge.species.taxonomy import TaxonomyRecord


def test_taxonomy_record():

    record = TaxonomyRecord(

        taxid="8128",

        scientific_name="Oreochromis niloticus",

        rank="species",

        lineage=["Eukaryota", "Chordata"],

    )

    assert record.taxid == "8128"

    assert record.rank == "species"

    assert record.scientific_name == "Oreochromis niloticus"

    assert "Chordata" in record.lineage
