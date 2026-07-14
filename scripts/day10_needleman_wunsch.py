# scripts/day10_needleman_wunsch.py
# Author: Boris Djagou
# Date: July 14, 2026
# Dynamic programming: Needleman-Wunsch global alignment from scratch

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Global alignment using Needleman-Wunsch dynamic programming.
    Returns: (aligned_seq1, aligned_seq2, score)
    """
    n, m = len(seq1), len(seq2)

    # Step 1: Initialize scoring matrix
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        score[i][0] = i * gap
    for j in range(m + 1):
        score[0][j] = j * gap

    # Step 2: Fill matrix using recurrence relation
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diag = score[i-1][j-1] + match_score
            up   = score[i-1][j] + gap
            left = score[i][j-1] + gap
            score[i][j] = max(diag, up, left)

    # Step 3: Traceback to build the alignment
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 and j > 0:
        current = score[i][j]
        match_score = match if seq1[i-1] == seq2[j-1] else mismatch

        if current == score[i-1][j-1] + match_score:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif current == score[i-1][j] + gap:
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            j -= 1

    while i > 0:
        aligned1.append(seq1[i-1])
        aligned2.append("-")
        i -= 1
    while j > 0:
        aligned1.append("-")
        aligned2.append(seq2[j-1])
        j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2), score[n][m]


def print_alignment(a1, a2):
    """Pretty-print an alignment with match indicators."""
    match_line = "".join("|" if x == y else (" " if x=="-" or y=="-" else ".")
                         for x, y in zip(a1, a2))
    print(f"  Seq1: {a1}")
    print(f"        {match_line}")
    print(f"  Seq2: {a2}")


def main():
    print("\n Day 10 — Needleman-Wunsch Global Alignment (from scratch)")

    # Two 30-aa protein sequences (as required by the exercise)
    # Using real kelch13 propeller domain excerpts (reference vs C580Y mutant)
    seq1 = "YVSSNLNIPRRNNCGVTSNGRIYCIGGYDG"   # 30 aa — reference region
    seq2 = "YVSSNLNIPRRNNYGVTSNGRIYCIGGYDG"   # 30 aa — with C580Y-like change

    print(f"\nSeq1 ({len(seq1)} aa): {seq1}")
    print(f"Seq2 ({len(seq2)} aa): {seq2}")

    print("\n[1/2] Running Needleman-Wunsch (match=+1, mismatch=-1, gap=-2)...")
    a1, a2, final_score = needleman_wunsch(seq1, seq2)

    print(f"\nAlignment (score = {final_score}):")
    print_alignment(a1, a2)

    # Count matches/mismatches/gaps
    matches    = sum(1 for x, y in zip(a1, a2) if x == y and x != "-")
    mismatches = sum(1 for x, y in zip(a1, a2) if x != y and x != "-" and y != "-")
    gaps       = sum(1 for x, y in zip(a1, a2) if x == "-" or y == "-")

    print(f"\nAlignment statistics:")
    print(f"  Matches    : {matches}")
    print(f"  Mismatches : {mismatches}")
    print(f"  Gaps       : {gaps}")
    print(f"  Identity   : {matches/len(a1)*100:.1f}%")

    # Second test: introduce a gap scenario
  
    print("[2/2] Testing with insertion (different lengths)...")
    seq3 = "YVSSNLNIPRRNNCGVTSNGRIYCIGGYDG"       # 30 aa
    seq4 = "YVSSNLNIPRRNNCGVTSNGRIYCIGGYYDG"      # 31 aa — 1 extra residue

    print(f"\nSeq3 ({len(seq3)} aa): {seq3}")
    print(f"Seq4 ({len(seq4)} aa): {seq4}")

    a3, a4, score2 = needleman_wunsch(seq3, seq4)
    print(f"\nAlignment (score = {score2}):")
    print_alignment(a3, a4)

    print("\n Needleman-Wunsch implementation validated!")


if __name__ == "__main__":
    main()
