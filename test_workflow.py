from pathlib import Path

from primerforge.config.config import Config
from primerforge.species.workflow import SpeciesWorkflow


def main():

    config = Config(
        Path("configs/example.yaml"),
    )

    workflow = SpeciesWorkflow(
        email=config.ncbi_email,
        api_key=config.ncbi_api_key,
    )

    result = workflow.run(
        species=config.organism,
        marker=config.gene,
        output_directory=Path("results"),
    )

    print()

    print("=" * 60)
    print("Workflow completed")
    print("=" * 60)

    print()

    print("Target FASTA:")
    print(result.dataset.target_fasta)

    print()

    print("Background FASTA:")
    print(result.dataset.background_fasta)

    print()

    print("Target alignment:")
    print(result.alignment.target_alignment)

    print()

    print("Background alignment:")
    print(result.alignment.background_alignment)

    print()

    print("Diagnostic SNPs:")
    print(len(result.diagnostic_sites))

    print()

    print("Diagnostic windows:")
    print(len(result.diagnostic_windows))

    if result.diagnostic_windows:

        window = result.diagnostic_windows[0]

        print()

        print("Best window")

        print(window.start)

        print(window.end)

        print(len(window.diagnostic_sites))


if __name__ == "__main__":

    main()
