"""
Tests for FASTAMerger.
"""

from pathlib import Path

from Bio import SeqIO

from primerforge.io.merge import FASTAMerger


def test_merge(tmp_path: Path):

    file1 = tmp_path / "a.fasta"

    file2 = tmp_path / "b.fasta"

    file1.write_text(
        ">a\nAAAA\n>b\nCCCC\n"
    )

    file2.write_text(
        ">c\nCCCC\n>d\nGGGG\n"
    )

    output = tmp_path / "merged.fasta"

    merger = FASTAMerger()

    merger.merge(
        [file1, file2],
        output,
    )

    records = list(
        SeqIO.parse(
            output,
            "fasta",
        )
    )

    assert len(records) == 3

    sequences = {
        str(r.seq)
        for r in records
    }

    assert sequences == {
        "AAAA",
        "CCCC",
        "GGGG",
    }
