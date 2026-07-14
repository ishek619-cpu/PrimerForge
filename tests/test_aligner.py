from pathlib import Path

from envoprimer.species.aligner import MAFFTAligner


def test_aligner():

    aligner = MAFFTAligner()

    assert aligner.executable == "mafft"
