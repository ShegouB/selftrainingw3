# scripts/day9_fetch_orthologs_v2.py
# Author: Boris Djagou
# Date: July 13, 2026
# Fixed: search-verified kelch13 orthologs instead of guessed accessions

from Bio import Entrez
import time
import os

Entrez.email = "djagouboris@gmail.com"

# Search by species + protein family name instead of hardcoded accessions
SPECIES_QUERIES = [
    {"species": "P. falciparum", "query": "kelch protein K13[Protein Name] AND Plasmodium falciparum[Organism]"},
    {"species": "P. vivax",      "query": "kelch[Protein Name] AND Plasmodium vivax[Organism] AND 600:800[SLEN]"},
    {"species": "P. knowlesi",   "query": "kelch[Protein Name] AND Plasmodium knowlesi[Organism] AND 600:800[SLEN]"},
    {"species": "P. malariae",   "query": "kelch[Protein Name] AND Plasmodium malariae[Organism] AND 600:800[SLEN]"},
    {"species": "P. ovale",      "query": "kelch[Protein Name] AND Plasmodium ovale[Organism] AND 600:800[SLEN]"},
]


def search_and_fetch(query, expected_min=600, expected_max=800):
    """Search NCBI protein DB, fetch first result matching expected length range."""
    handle = Entrez.esearch(db="protein", term=query, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]

    if not ids:
        return None, None, 0

    # Try each candidate until one matches expected length
    for pid in ids:
        handle = Entrez.efetch(db="protein", id=pid,
                               rettype="fasta", retmode="text")
        fasta = handle.read()
        handle.close()

        lines = fasta.strip().split("\n")
        seq   = "".join(lines[1:])
        length = len(seq)

        if expected_min <= length <= expected_max:
            header = lines[0]
            return header, seq, length

    # No match found within range — return first result anyway with warning
    handle = Entrez.efetch(db="protein", id=ids[0],
                           rettype="fasta", retmode="text")
    fasta = handle.read()
    handle.close()
    lines = fasta.strip().split("\n")
    seq   = "".join(lines[1:])
    return lines[0], seq, len(seq)


def main():
    print("\n Fixing kelch13 orthologs — search-verified accessions")

    verified_sequences = []

    for sp in SPECIES_QUERIES:
        print(f"\n[{sp['species']}]")
        header, seq, length = search_and_fetch(sp["query"])

        if seq is None:
            print(f"No results found")
            continue

        status = "OK" if 600 <= length <= 800 else "WARMING"
        print(f"  {status} {header[:70]}")
        print(f"  Length: {length} aa")

        verified_sequences.append({
            "species": sp["species"],
            "header":  header,
            "seq":     seq,
            "length":  length,
        })
        time.sleep(0.4)

    for v in verified_sequences:
        flag = "OKAY" if 600 <= v["length"] <= 800 else "CHECK"
        print(f"  {v['species']:<18} {v['length']:>5} aa  {flag}")

    # Save only the verified ones (within expected biological range)
    good = [v for v in verified_sequences if 600 <= v["length"] <= 800]
    print(f"\n{len(good)}/{len(verified_sequences)} orthologs within expected range (600-800 aa)")

    return good


if __name__ == "__main__":
    main()
