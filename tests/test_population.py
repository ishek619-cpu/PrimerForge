"""
Tests for population analysis.
"""

from pathlib import Path

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

from primerforge.analysis.population import PopulationAnalyzer


def test_population(tmp_path: Path):

    alignment = MultipleSeqAlignment([

        SeqRecord(
            Seq("ATGCGTACGT"),
            id="A",
        ),

        SeqRecord(
            Seq("ATGCGTACGT"),
            id="B",
        ),

        SeqRecord(
            Seq("ATGCGTACGT"),
            id="C",
        ),

    ])

    fasta = tmp_path / "alignment.fasta"

    AlignIO.write(
        alignment,
        fasta,
        "fasta",
    )

    analyzer = PopulationAnalyzer()

    result = analyzer.analyse(

        primer="ATGCG",

        alignment=fasta,

        start=0,

    )

    assert result["population_conservation"] == 100.0

    assert result["population_coverage"] == 100.0

    assert result["estimated_failure_rate"] == 0.0
