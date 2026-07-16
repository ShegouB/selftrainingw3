# scripts/day13_hmmscan.py
# Author: Boris Djagou
# Date: July 17, 2026
# Search the kelch13 profile HMM against a target sequence

import subprocess
import os


def run_hmmscan(hmm_db, query_fasta, output_file):
    """Run hmmscan to search a query against an HMM profile database."""
    cmd = ["hmmscan", "--tblout", output_file, hmm_db, query_fasta]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr


def run_hmmsearch(hmm_profile, target_db, output_file):
    """Run hmmsearch: search one HMM profile against a database of sequences."""
    cmd = ["hmmsearch", "--tblout", output_file, hmm_profile, target_db]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr


def main():
    print("\nDay 13 - HMM Search Against kelch13 Profile")
    print("=" * 60)

    hmm_profile = "results/hmm/kelch13.hmm"
    target_db = "data/kelch13_diverse.fasta"
    output_file = "results/hmm/hmmsearch_hits.tbl"

    if not os.path.exists(hmm_profile):
        print(f"HMM profile not found: {hmm_profile}")
        print("Run: hmmbuild results/hmm/kelch13.hmm results/alignments/kelch13_muscle.fasta")
        return

    print(f"\n[1/2] Searching profile HMM against target database...")
    print(f"  Profile: {hmm_profile}")
    print(f"  Targets: {target_db}")

    ok, out, err = run_hmmsearch(hmm_profile, target_db, output_file)

    if not ok:
        print(f"  hmmsearch failed: {err[:300]}")
        return

    print(f"  Results saved: {output_file}")

    print(f"\n[2/2] Parsing hit table...")
    print("-" * 60)

    if os.path.exists(output_file):
        with open(output_file) as f:
            lines = f.readlines()

        data_lines = [l for l in lines if not l.startswith("#")]

        if not data_lines:
            print("  No hits found.")
        else:
            print(f"  {'Target name':<30} {'E-value':>12} {'Score':>8}")
            print("  " + "-" * 52)
            for line in data_lines:
                fields = line.split()
                target_name = fields[0]
                evalue = fields[4]
                score = fields[5]
                print(f"  {target_name:<30} {evalue:>12} {score:>8}")

    print("\nHMM search complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
