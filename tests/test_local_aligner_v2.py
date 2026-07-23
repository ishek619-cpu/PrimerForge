"""
Unit test for EnvoPrimer v2 LocalAligner.
"""

from envoprimer_v2.pcr.local_aligner import LocalAligner


def main():

    primer = "TGGAGGCTTTGGAAACTGAC"

    window = (
        "AAA"
        "TGGAGGCTTTGGAAACTGAC"
        "TTTTTTTTTT"
    )

    aligner = LocalAligner()

    result = aligner.align(

        primer=primer,

        window=window,

        sequence_id=0,

        strand="+",

        window_start=0,

    )

    print()

    print("=" * 60)
    print("LOCAL ALIGNMENT TEST")
    print("=" * 60)

    print("Start :", result.start)
    print("End   :", result.end)

    print("Identity :", result.identity)
    print("Mismatch :", result.mismatches)
    print("Terminal :", result.terminal_mismatches)

    print("Score :", result.score)

    print()

    print(result.primer_alignment)
    print(result.template_alignment)


if __name__ == "__main__":
    main()
