# scripts/day9_fetch_kelch13_strains.py
# Author: Boris Djagou
# Date: July 13, 2026
# Fetch kelch13 protein sequences from multiple P. falciparum strains

from Bio import Entrez, SeqIO
import time
import os

Entrez.email = "djagouboris@gmail.com"

# Search for kelch13 sequences from different geographic origins
QUERIES = [
    ("kelch13_diverse",
     "kelch13[gene] AND Plasmodium falciparum[organism] AND 600:800[slen]",
     15),
]

def fetch_sequences(query, max_results=15):
    """Search and fetch protein sequences."""
    handle = Entrez.esearch(db="protein", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    print(f"  Found {len(ids)} sequences")

    if not ids:
        return []

    handle = Entrez.efetch(db="protein", id=",".join(ids),
                           rettype="fasta", retmode="text")
    fasta_text = handle.read()
    handle.close()
    return fasta_text


def main():
    os.makedirs("data", exist_ok=True)
    print("\n Day 9 — Fetching kelch13 sequences from multiple strains")


    for name, query, n in QUERIES:
        print(f"\nQuery: {query[:60]}...")
        fasta = fetch_sequences(query, max_results=n)

        if fasta:
            path = f"data/{name}.fasta"
            with open(path, "w") as f:
                f.write(fasta)

            # Count and show sequences
            seqs = [l for l in fasta.split("\n") if l.startswith(">")]
            print(f"  Sequences retrieved: {len(seqs)}")
            print(f"  Saved: {path}")
            print("\n  Sequence headers:")
            for s in seqs[:10]:
                print(f"    {s[:75]}")
        time.sleep(0.5)

    print("\n Sequences ready for alignment")

if __name__ == "__main__":
    main()
