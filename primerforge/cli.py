"""
PrimerForge command line interface.
"""

from pathlib import Path

import typer

from primerforge.core.pipeline import Pipeline

app = typer.Typer(
    add_completion=False,
)


@app.command()
def run(
    config: Path,
):
    """
    Run PrimerForge.
    """

    pipeline = Pipeline()

    pipeline.run(config)


def main():

    app()


if __name__ == "__main__":

    main()
