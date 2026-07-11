from pathlib import Path

from primerforge.species.diagnostics import DiagnosticFinder


def test_finder(tmp_path: Path):

    target = tmp_path / "target.fasta"

    background = tmp_path / "background.fasta"

    target.write_text(
""">a
AAAA
>b
AAAA
""",
encoding="utf-8",
)

    background.write_text(
""">c
AATA
>d
AATA
""",
encoding="utf-8",
)

    finder = DiagnosticFinder()

    sites = finder.find(
        target,
        background,
    )

    assert len(sites) == 1

    assert sites[0].position == 2
