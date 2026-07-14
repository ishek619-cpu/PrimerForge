"""
Primer3 interface for EnvoPrimer.

Designs primers from a diagnostic window.
"""

from __future__ import annotations

import primer3

from envoprimer_v2.windows.builder import DiagnosticWindow


class PrimerDesigner:

    def __init__(

        self,

        optimum_tm: float = 60.0,

        minimum_tm: float = 58.0,

        maximum_tm: float = 62.0,

    ):

        self.optimum_tm = optimum_tm
        self.minimum_tm = minimum_tm
        self.maximum_tm = maximum_tm

    def design(

        self,

        window: DiagnosticWindow,

    ):

        sequence_arguments = {

            "SEQUENCE_TEMPLATE": window.sequence,

        }

        global_arguments = {

            #
            # Primer length
            #
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 25,

            #
            # Tm
            #
            "PRIMER_OPT_TM": self.optimum_tm,
            "PRIMER_MIN_TM": self.minimum_tm,
            "PRIMER_MAX_TM": self.maximum_tm,

            #
            # GC
            #
            "PRIMER_MIN_GC": 40.0,
            "PRIMER_MAX_GC": 60.0,

            #
            # eDNA amplicon
            #
            "PRIMER_PRODUCT_SIZE_RANGE": [[70, 200]],

            #
            # Return many candidates
            #
            "PRIMER_NUM_RETURN": 100,

            #
            # Basic quality filters
            #
            "PRIMER_MAX_POLY_X": 4,

            "PRIMER_SALT_MONOVALENT": 50.0,

            "PRIMER_DNA_CONC": 50.0,

        }

        return primer3.bindings.design_primers(

            seq_args=sequence_arguments,

            global_args=global_arguments,

        )
