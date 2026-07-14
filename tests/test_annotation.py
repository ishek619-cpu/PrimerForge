"""
Tests for GenBank annotation parsing.
"""

from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import (
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from envoprimer.reference.annotation import (
    AnnotationParser,
)


def create_genbank(path: Path):

    record = SeqRecord(
        Seq("ATGC" * 500),
        id="NC_TEST",
        name="NC_TEST",
        description="Synthetic mitochondrial genome",
    )

    #
    # Required for writing GenBank files.
    #
    record.annotations["molecule_type"] = "DNA"

    record.features = [

        SeqFeature(
            FeatureLocation(
                0,
                100,
                strand=1,
            ),
            type="gene",
            qualifiers={
                "gene": ["COI"],
            },
        ),

        SeqFeature(
            FeatureLocation(
                150,
                350,
                strand=1,
            ),
            type="CDS",
            qualifiers={
                "gene": ["COI"],
                "product": ["cytochrome oxidase I"],
            },
        ),

        SeqFeature(
            FeatureLocation(
                500,
                650,
                strand=-1,
            ),
            type="gene",
            qualifiers={
                "gene": ["CYTB"],
            },
        ),

        SeqFeature(
            FeatureLocation(
                700,
                850,
                strand=1,
            ),
            type="rRNA",
            qualifiers={
                "gene": ["12S"],
            },
        ),

    ]

    SeqIO.write(
        record,
        path,
        "genbank",
    )


def test_load_annotation(tmp_path: Path):

    gb = tmp_path / "test.gb"

    create_genbank(
        gb,
    )

    parser = AnnotationParser()

    genes = parser.load(
        gb,
    )

    assert len(genes) == 4


def test_find_gene(tmp_path: Path):

    gb = tmp_path / "test.gb"

    create_genbank(
        gb,
    )

    parser = AnnotationParser()

    parser.load(
        gb,
    )

    gene = parser.get_gene(
        "COI",
    )

    assert gene is not None

    assert gene.name == "COI"

    assert gene.start == 1

    assert gene.end == 100


def test_missing_gene(tmp_path: Path):

    gb = tmp_path / "test.gb"

    create_genbank(
        gb,
    )

    parser = AnnotationParser()

    parser.load(
        gb,
    )

    assert parser.get_gene(
        "ND6",
    ) is None


def test_feature_types(tmp_path: Path):

    gb = tmp_path / "test.gb"

    create_genbank(
        gb,
    )

    parser = AnnotationParser()

    parser.load(
        gb,
    )

    kinds = parser.feature_types()

    assert "gene" in kinds

    assert "CDS" in kinds

    assert "rRNA" in kinds


def test_gene_names(tmp_path: Path):

    gb = tmp_path / "test.gb"

    create_genbank(
        gb,
    )

    parser = AnnotationParser()

    parser.load(
        gb,
    )

    names = parser.gene_names()

    assert "COI" in names

    assert "CYTB" in names

    assert "12S" in names


def test_gene_length():

    from envoprimer.reference.annotation import Gene

    gene = Gene(

        name="COI",

        feature_type="gene",

        start=100,

        end=500,

        strand=1,

    )

    assert gene.length == 401


def test_gene_contains():

    from envoprimer.reference.annotation import Gene

    gene = Gene(

        name="COI",

        feature_type="gene",

        start=100,

        end=500,

        strand=1,

    )

    assert gene.contains(100)

    assert gene.contains(250)

    assert gene.contains(500)

    assert not gene.contains(99)

    assert not gene.contains(501)


def test_gene_to_dict():

    from envoprimer.reference.annotation import Gene

    gene = Gene(

        name="COI",

        feature_type="gene",

        start=10,

        end=20,

        strand=1,

    )

    data = gene.to_dict()

    assert data["name"] == "COI"

    assert data["type"] == "gene"

    assert data["length"] == 11
