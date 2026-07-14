"""
EnvoPrimer command line interface.
"""

from __future__ import annotations

import argparse


def build_parser():

    parser = argparse.ArgumentParser(

        prog="envoprimer",

        description="EnvoPrimer species-specific eDNA primer designer.",

    )

    parser.add_argument(

        "--config",

        required=True,

        help="Configuration YAML file.",

    )

    return parser


def main():

    parser = build_parser()

    args = parser.parse_args()

    print()

    print("=" * 70)
    print("EnvoPrimer")
    print("=" * 70)

    print()

    print("Configuration")

    print(args.config)

    print()

    print("Pipeline integration begins in next step.")

    return 0


if __name__ == "__main__":

    raise SystemExit(main())
