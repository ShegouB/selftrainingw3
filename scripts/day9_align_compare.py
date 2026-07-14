# scripts/day9_align_compare.py
# Author: Boris Djagou
# Date: July 13, 2026
# Fixed: MUSCLE v5 syntax, proper subprocess list args, robust error handling

import subprocess
import os
import time
from Bio import AlignIO

INPUT   = "data/kelch13_diverse.fasta"
OUT_DIR = "results/alignments"
os.makedirs(OUT_DIR, exist_ok=True)


def check_input():
    """Verify input file exists and has sequences before running aligners."""
    if not os.path.exists(INPUT):
        print(f"Input file not found: {INPUT}")
        print("Run day9_fetch_kelch13_strains.py first.")
        return False
    with open(INPUT) as f:
        content = f.read()
    n_seqs = content.count(">")
    if n_seqs < 2:
        print(f" Only {n_seqs} sequence(s) found — need at least 2 to align.")
        return False
    print(f"Input verified: {n_seqs} sequences in {INPUT}")
    return True


def run_muscle(input_fasta, output_afa):
    """Run MUSCLE v5 alignment — uses -align/-output syntax."""
    cmd = ["muscle", "-align", input_fasta, "-output", output_afa]
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start
    if result.returncode == 0 and os.path.exists(output_afa):
        print(f"  MUSCLE   : {elapsed:.2f}s → {output_afa}")
        return True
    print(f"  MUSCLE   : {result.stderr[:300]}")
    return False


def run_mafft(input_fasta, output_afa):
    """Run MAFFT alignment — writes to stdout, redirect to file."""
    cmd = ["mafft", "--auto", "--quiet", input_fasta]
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start
    if result.returncode == 0 and result.stdout.strip():
        with open(output_afa, "w") as f:
            f.write(result.stdout)
        print(f"  MAFFT    : {elapsed:.2f}s → {output_afa}")
        return True
    print(f"  MAFFT    : {result.stderr[:300]}")
    return False


def run_clustalw(input_fasta, output_aln):
    """Run ClustalW alignment."""
    cmd = ["clustalw2", f"-INFILE={input_fasta}",
           f"-OUTFILE={output_aln}", "-OUTPUT=FASTA", "-QUIET"]
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start
    if result.returncode == 0 and os.path.exists(output_aln):
        print(f"  ClustalW : {elapsed:.2f}s → {output_aln}")
        return True
    print(f"  ClustalW : {result.stderr[:300]}")
    return False


def analyse_alignment(filepath, name, fmt="fasta"):
    """Analyse an alignment: length, gaps, conservation."""
    if not os.path.exists(filepath):
        print(f"\n  [{name}] — file not found, skipping analysis")
        return
    try:
        aln = AlignIO.read(filepath, fmt)
    except Exception as e:
        print(f"\n  [{name}] — could not read alignment ({e})")
        return

    n_seq   = len(aln)
    aln_len = aln.get_alignment_length()

    conserved = 0
    gap_cols  = 0
    for i in range(aln_len):
        col = aln[:, i]
        if "-" in col:
            gap_cols += 1
        if len(set(col)) == 1 and "-" not in col:
            conserved += 1

    pct_conserved = conserved / aln_len * 100
    pct_gaps      = gap_cols / aln_len * 100

    print(f"\n  [{name}]")
    print(f"    Sequences      : {n_seq}")
    print(f"    Alignment len  : {aln_len:,} positions")
    print(f"    Conserved cols : {conserved:,} ({pct_conserved:.1f}%)")
    print(f"    Gap columns    : {gap_cols:,} ({pct_gaps:.1f}%)")

    # Check conservation around known resistance domain (kelch propeller ~440-726)
    start_pos = min(440, aln_len - 1)
    end_pos   = min(680, aln_len)
    resist_conserved = 0
    for i in range(start_pos, end_pos):
        col = aln[:, i]
        if len(set(col)) == 1 and "-" not in col:
            resist_conserved += 1
    print(f"    Conserved in propeller domain (440-680): {resist_conserved} positions")


def main():
    print("\n Day 9 — Multiple Sequence Alignment Comparison")

    if not check_input():
        return

    muscle_out  = f"{OUT_DIR}/kelch13_muscle.fasta"
    mafft_out   = f"{OUT_DIR}/kelch13_mafft.fasta"
    clustal_out = f"{OUT_DIR}/kelch13_clustalw.fasta"

    print("\n[1/2] Running aligners...")
    run_muscle(INPUT, muscle_out)
    run_mafft(INPUT, mafft_out)
    run_clustalw(INPUT, clustal_out)

    print("\n[2/2] Analysing alignments...")
    for path, name in [
        (muscle_out,  "MUSCLE"),
        (mafft_out,   "MAFFT"),
        (clustal_out, "ClustalW"),
    ]:
        analyse_alignment(path, name)

    print("\n Alignment comparison complete!")
    print(f"   Output files in: {OUT_DIR}/")

if __name__ == "__main__":
    main()
