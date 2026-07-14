# scripts/day9_build_dataset_v3.py
# Author: Boris Djagou
# Date: July 13, 2026
# Final version: verified P. falciparum kelch13 reference + 9 WHO resistance mutants

from Bio import Entrez
import time
import os

Entrez.email = "djagouboris@gmail.com"

# Verified accession from Week 2 — confirmed 726 aa
REFERENCE_ACCESSION = "XP_001350158.1"

VALIDATED_MUTATIONS = [
    {"name": "F446I",  "pos": 446, "wt": "F", "mut": "I", "region": "Myanmar - highly prevalent"},
    {"name": "N458Y",  "pos": 458, "wt": "N", "mut": "Y", "region": "Southeast Asia"},
    {"name": "M476I",  "pos": 476, "wt": "M", "mut": "I", "region": "Southeast Asia"},
    {"name": "Y493H",  "pos": 493, "wt": "Y", "mut": "H", "region": "Southeast Asia"},
    {"name": "R539T",  "pos": 539, "wt": "R", "mut": "T", "region": "Southeast Asia, emerging in Africa"},
    {"name": "I543T",  "pos": 543, "wt": "I", "mut": "T", "region": "Southeast Asia"},
    {"name": "P553L",  "pos": 553, "wt": "P", "mut": "L", "region": "Southeast Asia"},
    {"name": "R561H",  "pos": 561, "wt": "R", "mut": "H", "region": "Rwanda, DR Congo — first African emergence"},
    {"name": "C580Y",  "pos": 580, "wt": "C", "mut": "Y", "region": "Global — most common resistance mutation"},
]


def fetch_reference():
    """Fetch the verified P. falciparum kelch13 reference sequence."""
    handle = Entrez.efetch(db="protein", id=REFERENCE_ACCESSION,
                           rettype="fasta", retmode="text")
    fasta = handle.read()
    handle.close()
    lines  = fasta.strip().split("\n")
    seq    = "".join(lines[1:])
    return seq


def build_mutant(reference_seq, mutation):
    """Apply a validated point mutation to the reference sequence."""
    pos = mutation["pos"] - 1  # 0-indexed
    actual_wt = reference_seq[pos]

    match = "OK" if actual_wt == mutation["wt"] else "WARMINGS"
    print(f"  {mutation['name']:<8} pos {mutation['pos']:>4}  "
          f"expected {mutation['wt']}, found {actual_wt}  {match}")

    mutant = reference_seq[:pos] + mutation["mut"] + reference_seq[pos+1:]
    return mutant


def main():
    os.makedirs("data", exist_ok=True)
    print("\n Day 9 v3 — Verified kelch13 Dataset: Reference + 9 WHO Mutants")

    print(f"\nFetching reference: {REFERENCE_ACCESSION}")
    reference_seq = fetch_reference()
    print(f"Reference length: {len(reference_seq)} aa")

    all_fasta = [f">Pf_kelch13_REFERENCE_3D7\n{reference_seq}"]

    print(f"\nBuilding {len(VALIDATED_MUTATIONS)} WHO-validated resistance mutants:")

    for mut in VALIDATED_MUTATIONS:
        mutant_seq = build_mutant(reference_seq, mut)
        header = f">Pf_kelch13_{mut['name']}"
        all_fasta.append(f"{header}\n{mutant_seq}")

    output_path = "data/kelch13_diverse.fasta"
    with open(output_path, "w") as f:
        f.write("\n".join(all_fasta) + "\n")

    print(f"\n Dataset built: {len(all_fasta)} sequences "
          f"(1 reference + {len(VALIDATED_MUTATIONS)} mutants)")
    print(f"   All sequences: {len(reference_seq)} aa (identical length except mutation point)")
    print(f"   Saved: {output_path}")


if __name__ == "__main__":
    main()
