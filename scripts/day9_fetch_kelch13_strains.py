# scripts/day9_fetch_kelch13_strains.py
# Author: Boris Djagou
# Date: July 13, 2026
# Fixed: kelch13 orthologs across Plasmodium species + WHO-validated resistance mutants

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import time
import os

Entrez.email = "djagouboris@gmail.com"

# ── PART A: kelch13 orthologs across Plasmodium species (real NCBI data) ────
ORTHOLOGS = [
    {"species": "P. falciparum", "accession": "XP_001350158.1"},
    {"species": "P. vivax",      "accession": "XP_001615711.1"},
    {"species": "P. knowlesi",   "accession": "XP_002259541.1"},
    {"species": "P. malariae",   "accession": "SBT77415.1"},
    {"species": "P. reichenowi", "accession": "XP_002805349.1"},
]

# ── PART B: WHO-validated kelch13 resistance mutations (literature-verified) ─
# Source: WHO artemisinin resistance mutation catalogue
# 9 validated mutations that confer confirmed artemisinin resistance
VALIDATED_MUTATIONS = [
    {"name": "F446I",  "pos": 446, "wt": "F", "mut": "I", "region": "Myanmar - highly prevalent"},
    {"name": "N458Y",  "pos": 458, "wt": "N", "mut": "Y", "region": "Southeast Asia"},
    {"name": "M476I",  "pos": 476, "wt": "M", "mut": "I", "region": "Southeast Asia"},
    {"name": "Y493H",  "pos": 493, "wt": "Y", "mut": "H", "region": "Southeast Asia"},
    {"name": "R539T",  "pos": 539, "wt": "R", "mut": "T", "region": "Southeast Asia, emerging in Africa"},
    {"name": "I543T",  "pos": 543, "wt": "I", "mut": "T", "region": "Southeast Asia"},
    {"name": "P553L",  "pos": 553, "wt": "P", "mut": "L", "region": "Southeast Asia"},
    {"name": "C580Y",  "pos": 580, "wt": "C", "mut": "Y", "region": "Global — most common resistance mutation"},
    {"name": "R561H",  "pos": 561, "wt": "R", "mut": "H", "region": "Rwanda, DR Congo — first African emergence"},
]


def fetch_ortholog(species, accession):
    """Fetch a kelch13 ortholog protein sequence."""
    try:
        handle = Entrez.efetch(db="protein", id=accession,
                               rettype="fasta", retmode="text")
        fasta = handle.read()
        handle.close()
        return fasta
    except Exception as e:
        print(f"    ❌ Failed to fetch {accession}: {e}")
        return None


def build_mutant_sequence(reference_seq, mutation):
    """Apply a single point mutation to the reference sequence."""
    pos = mutation["pos"] - 1  # 0-indexed
    if pos >= len(reference_seq):
        return None
    if reference_seq[pos] != mutation["wt"]:
        print(f"    ⚠ Position {mutation['pos']} is "
              f"{reference_seq[pos]} not {mutation['wt']} — check numbering")
    mutant = reference_seq[:pos] + mutation["mut"] + reference_seq[pos+1:]
    return mutant


def main():
    os.makedirs("data", exist_ok=True)
    print("\n🧬 Day 9 — kelch13 Diversity Dataset (Orthologs + Resistance Mutants)")
    print("=" * 70)

    all_fasta = []

    # ── PART A: Fetch orthologs ──────────────────────────────────────
    print("\n[1/2] Fetching kelch13 orthologs across Plasmodium species")
    print("-" * 70)

    reference_seq = None
    for orth in ORTHOLOGS:
        print(f"\n  {orth['species']} ({orth['accession']})")
        fasta = fetch_ortholog(orth["species"], orth["accession"])

        if fasta and len(fasta) > 20:
            lines  = fasta.strip().split("\n")
            header = f">{orth['species'].replace(' ', '_')}_{orth['accession']}"
            seq    = "".join(lines[1:])
            print(f"    ✅ {len(seq)} aa")
            all_fasta.append(f"{header}\n{seq}")

            if orth["species"] == "P. falciparum":
                reference_seq = seq
        else:
            print(f"    ❌ Not retrieved")
        time.sleep(0.4)

    # ── PART B: Build WHO-validated resistance mutants ──────────────
    print(f"\n\n[2/2] Building 9 WHO-validated resistance mutants")
    print("-" * 70)

    if reference_seq:
        for mut in VALIDATED_MUTATIONS:
            mutant_seq = build_mutant_sequence(reference_seq, mut)
            if mutant_seq:
                header = f">Pf_kelch13_{mut['name']}_mutant"
                all_fasta.append(f"{header}\n{mutant_seq}")
                print(f"  {mut['name']:<8} pos {mut['pos']:>4}  "
                      f"{mut['wt']}→{mut['mut']}   {mut['region']}")
    else:
        print("  ❌ No P. falciparum reference sequence available for mutagenesis")

    # ── Save combined FASTA ───────────────────────────────────────────
    output_path = "data/kelch13_diverse.fasta"
    with open(output_path, "w") as f:
        f.write("\n".join(all_fasta) + "\n")

    n_seqs = len(all_fasta)
    print(f"\n\n✅ Dataset built: {n_seqs} sequences")
    print(f"   ({len(ORTHOLOGS)} orthologs + {len(VALIDATED_MUTATIONS)} resistance mutants)")
    print(f"   Saved: {output_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
