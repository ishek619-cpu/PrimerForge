"""
Tests for PrimerForge configuration.
"""

from pathlib import Path

from primerforge.config.config import Config


def test_configuration():

    cfg = Config(
        Path(
            "configs/example.yaml"
        )
    )

    assert cfg.project != ""

    assert cfg.organism == "Oreochromis niloticus"

    assert cfg.gene == "CYTB"

    assert cfg.taxonomic_rank == "species"

    assert cfg.version == "0.1.0"


def test_dictionary_access():

    cfg = Config(
        Path(
            "configs/example.yaml"
        )
    )

    assert cfg["marker"]["gene"] == "CYTB"


def test_default_value():

    cfg = Config(
        Path(
            "configs/example.yaml"
        )
    )

    assert cfg.get(
        "xyz",
        "default",
    ) == "default"
