# scripts/day12_real_blast.py
# Author: Boris Djagou
# Date: July 16, 2026
# Run real BLAST+ locally and compare against simplified implementation

import subprocess
import os

DB_DIR = "results/blast_db"
os.makedirs(DB_DIR, exist_ok=True)


def build_blast_db(fasta_path, db_name):
    """Build a local BLAST protein database."""
    cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "prot",
           "-out", db_name, "-title", "kelch13_mutants"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr


def run_blastp(query_fasta, db_name, output_file, evalue=10):
    """Run blastp search against a local database."""
    cmd = ["blastp", "-query", query_fasta, "-db", db_name,
           "-out", output_file, "-outfmt",
           "6 qseqid sseqid pident length mismatch gapopen "
           "qstart qend sstart send evalue bitscore",
           "-evalue", str(evalue)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr


def main():
    print("\nDay 12 - Real BLASTP Search")

    fasta_db = "data/kelch13_diverse.fasta"
    db_name = f"{DB_DIR}/kelch13_db"

    print(f"\n[1/3] Building local BLAST database from {fasta_db}...")
    ok, out, err = build_blast_db(fasta_db, db_name)
    if not ok:
        print(f"  Database build failed: {err[:300]}")
        return
    print(f"  Database built: {db_name}")

    query_path = "data/query_mutated.fasta"
    with open(query_path, "w") as f:
        f.write(">mutated_query\nCIGGYDGSSIIPKVEAYAHRMKAWVEVAPL\n")
    print(f"\n[2/3] Query saved: {query_path}")

    output_path = "results/day12_blastp_hits.tsv"
    print(f"\n[3/3] Running blastp search (E-value threshold = 10)...")
    ok, out, err = run_blastp(query_path, db_name, output_path, evalue=10)

    if not ok:
        print(f"  BLASTP failed: {err[:300]}")
        return

    print(f"  Results saved: {output_path}")
    print("\nHit table (qseqid, sseqid, pident, length, evalue, bitscore):")

    if os.path.exists(output_path):
        with open(output_path) as f:
            lines = f.readlines()

        if not lines:
            print("  No hits found within E-value threshold.")
        else:
            print(f"  {'Subject':<20} {'%ident':>7} {'length':>7} "
                  f"{'evalue':>12} {'bitscore':>9}")
            for line in lines:
                fields = line.strip().split("\t")
                qseqid, sseqid, pident, length = fields[0:4]
                evalue, bitscore = fields[10], fields[11]
                print(f"  {sseqid:<20} {pident:>7} {length:>7} "
                      f"{evalue:>12} {bitscore:>9}")

    print("\nReal BLASTP search complete")


if __name__ == "__main__":
    main()
