"""
PrimerForge Command Line Interface
"""

from pathlib import Path

import typer

from primerforge import __version__
from primerforge.core.pipeline import Pipeline

app = typer.Typer(
    help="PrimerForge - Species-specific PCR primer discovery platform"
)


@app.command()
def version():
    """
    Show PrimerForge version.
    """
    typer.echo("PrimerForge")
    typer.echo(f"Version: {__version__}")


@app.command()
def run(config: Path):
    """
    Run the complete PrimerForge pipeline.
    """
    pipeline = Pipeline(config)
    pipeline.run()


def main():
    app()


if __name__ == "__main__":
    main()
