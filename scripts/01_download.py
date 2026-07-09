#!/usr/bin/env python3

"""
PrimerForge
Module 1 : Sequence Downloader

Downloads all GenBank records for a target species and marker.

Example:
python 01_download.py --target "Oreochromis niloticus" --marker CYTB
"""

import argparse
import os
from Bio import Entrez
from Bio import SeqIO

# -----------------------------

Entrez.email = "ishek619@gmail.com"

# -----------------------------

parser = argparse.ArgumentParser(
    description="Download sequences from NCBI"
)

parser.add_argument(
    "--target",
    required=True,
    help="Scientific name"
)

parser.add_argument(
    "--marker",
    required=True,
    help="Gene name (e.g. CYTB, COI, 12S)"
)

args = parser.parse_args()

species = args.target
marker = args.marker

print("=" * 60)
print("PrimerForge")
print("=" * 60)
print("Target :", species)
print("Marker :", marker)
print()

folder = os.path.join(
    "..",
    "data",
    "genes",
    species.replace(" ", "_")
)

os.makedirs(folder, exist_ok=True)

query = f'"{species}"[Organism] AND {marker}[Gene]'

print("Searching NCBI...")
print(query)

search = Entrez.esearch(
    db="nucleotide",
    term=query,
    retmax=10000
)

ids = Entrez.read(search)["IdList"]

print()
print(f"Found {len(ids)} records")

if len(ids) == 0:
    quit()

gb_file = os.path.join(folder, f"{marker}.gb")
fa_file = os.path.join(folder, f"{marker}.fasta")

records = []

for i, accession in enumerate(ids, 1):

    print(f"[{i}/{len(ids)}] {accession}")

    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    try:
        record = SeqIO.read(handle, "genbank")
        records.append(record)

    except:
        pass

SeqIO.write(records, gb_file, "genbank")
SeqIO.write(records, fa_file, "fasta")

print()
print("=" * 60)
print("Download complete")
print(gb_file)
print(fa_file)
print("=" * 60)
