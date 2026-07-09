"""
PrimerForge CLI

Entry point for all PrimerForge commands.
"""

from typing import Optional

import typer
from rich import print

app = typer.Typer(
    help="PrimerForge - Automated species-specific primer discovery"
)


@app.command()
def version():
    """
    Show PrimerForge version.
    """
    print("[bold green]PrimerForge[/bold green]")
    print("Version: 0.1.0")


@app.command()
def download(
    species: str = typer.Option(..., help="Target species"),
    gene: str = typer.Option(..., help="Target gene"),
):
    """
    Download sequences from NCBI.
    """
    print(f"[cyan]Species:[/cyan] {species}")
    print(f"[cyan]Gene:[/cyan] {gene}")
    print()
    print("Download module not implemented yet.")


if __name__ == "__main__":
    app()
