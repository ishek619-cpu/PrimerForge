"""
Biological species name filter.

Removes hybrids, unidentified taxa, strains,
environmental samples and other names unsuitable
for species-specific primer design.
"""

from __future__ import annotations


class SpeciesFilter:
    """
    Keep only valid binomial species names.
    """

    INVALID_TERMS = {

        " sp.",
        " cf.",
        " aff.",
        " hybrid",
        " strain",
        " clone",
        " isolate",
        " voucher",
        " uncultured",
        " environmental",
        " metagenome",
        " x ",
        " x",
        "x ",
        "BOLD",

    }

    def valid(
        self,
        species: str,
    ) -> bool:

        name = species.strip()

        #
        # Reject names containing unwanted terms.
        #
        lower = name.lower()

        for term in self.INVALID_TERMS:

            if term.lower() in lower:

                return False

        #
        # Must be exactly Genus species.
        #
        parts = name.split()

        if len(parts) != 2:

            return False

        #
        # Genus capitalized.
        #
        if not parts[0][0].isupper():

            return False

        #
        # Species epithet lowercase.
        #
        if not parts[1].islower():

            return False

        return True

    def filter(
        self,
        species_list: list[str],
    ) -> list[str]:

        kept = []

        seen = set()

        for species in species_list:

            if not self.valid(
                species,
            ):
                continue

            if species in seen:
                continue

            kept.append(
                species,
            )

            seen.add(
                species,
            )

        return sorted(
            kept,
        )
