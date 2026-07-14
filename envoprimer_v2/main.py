"""
EnvoPrimer executable entry point.
"""

from __future__ import annotations

from pathlib import Path

from envoprimer_v2.pipeline.workflow import EnvoPrimerWorkflow


def run(

    reference_sequence: str,

    target_alignment: Path,

    background_alignment: Path,

    background_sequences: dict,

    output_directory: Path = Path("results"),

):

    workflow = EnvoPrimerWorkflow()

    return workflow.run(

        reference_sequence=reference_sequence,

        target_alignment=target_alignment,

        background_alignment=background_alignment,

        background_sequences=background_sequences,

        output_directory=output_directory,

    )
