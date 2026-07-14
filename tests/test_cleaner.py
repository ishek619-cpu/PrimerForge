from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from envoprimer.species.cleaner import SequenceCleaner


def test_cleaner(tmp_path: Path):

    records = [

        SeqRecord(
            Seq("ATGC" * 30),
            id="1",
        ),

        SeqRecord(
            Seq("ATGC" * 30),
            id="2",
        ),

        SeqRecord(
            Seq("NNNN"),
            id="3",
        ),

    ]

    input_fasta = tmp_path / "input.fasta"

    output_fasta = tmp_path / "output.fasta"

    SeqIO.write(
        records,
        input_fasta,
        "fasta",
    )

    cleaner = SequenceCleaner()

    kept = cleaner.clean(
        input_fasta,
        output_fasta,
    )

    assert kept == 1
