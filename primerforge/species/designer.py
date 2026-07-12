"""
Species-specific primer designer.

Connects SpeciesWorkflow to PrimerDiscovery.
"""

from __future__ import annotations

from pathlib import Path

from primerforge.primer.discovery import PrimerDiscovery
from primerforge.specificity.engine import SpecificityEngine
from primerforge.species.workflow import WorkflowResult


class SpeciesPrimerDesigner:
    """
    Design primers from a SpeciesWorkflow result.
    """

    def __init__(self):

        self.discovery = PrimerDiscovery()

    def design(
        self,
        workflow: WorkflowResult,
        specificity_engine: SpecificityEngine | None = None,
    ):

        if not workflow.diagnostic_windows:

            raise RuntimeError(
                "No diagnostic windows were found."
            )

        #
        # Highest-ranked diagnostic window
        #
        window = workflow.diagnostic_windows[0]

        #
        # Convert to Primer3 region
        #
        region = window.to_region()

        #
        # Design primers
        #
        return self.discovery.discover(

            reference_fasta=workflow.dataset.target_fasta,

            alignment_fasta=workflow.alignment.target_alignment,

            specificity_engine=specificity_engine,

            region=region,

        )
