from pathlib import Path

from primerforge.species.aligner import MAFFTAligner


def test_aligner():

    aligner = MAFFTAligner()

    assert aligner.executable == "mafft"
