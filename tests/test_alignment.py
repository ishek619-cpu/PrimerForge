"""
Tests for MAFFT alignment.
"""

from pathlib import Path

from primerforge.analysis.align import MAFFTAligner


def test_alignment(tmp_path: Path):

    fasta = tmp_path / "input.fasta"

    fasta.write_text(
        ">a\nATGCTAGCTAGCTAGC\n"
        ">b\nATGCTAGCTAGCTTGC\n"
    )

    output = tmp_path / "alignment.fasta"

    aligner = MAFFTAligner()

    aligner.align(
        fasta,
        output,
    )

    assert output.exists()

    text = output.read_text()

    assert text.startswith(">")

    assert text.count(">") == 2
