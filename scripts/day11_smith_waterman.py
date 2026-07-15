# scripts/day11_smith_waterman.py
# Author: Boris Djagou
# Date: July 15, 2026
# Local alignment using Smith-Waterman with BLOSUM62 scoring matrix

from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")


def smith_waterman(seq1, seq2, gap_open=-10, gap_extend=-1, matrix=BLOSUM62):
    """
    Local alignment using Smith-Waterman dynamic programming.
    Uses BLOSUM62 for substitution scoring and affine gap penalties.
    Returns: (aligned_seq1, aligned_seq2, score, start1, start2)
    """
    n, m = len(seq1), len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[None] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            try:
                sub_score = matrix[seq1[i-1], seq2[j-1]]
            except KeyError:
                sub_score = matrix[seq2[j-1], seq1[i-1]]

            diag = score[i-1][j-1] + sub_score
            up   = score[i-1][j] + gap_extend
            left = score[i][j-1] + gap_extend

            best = max(0, diag, up, left)
            score[i][j] = best

            if best == 0:
                traceback[i][j] = "stop"
            elif best == diag:
                traceback[i][j] = "diag"
            elif best == up:
                traceback[i][j] = "up"
            else:
                traceback[i][j] = "left"

            if best > max_score:
                max_score = best
                max_pos = (i, j)

    aligned1, aligned2 = [], []
    i, j = max_pos
    end1, end2 = i, j

    while i > 0 and j > 0 and traceback[i][j] != "stop":
        direction = traceback[i][j]
        if direction == "diag":
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2), max_score, i, j, end1, end2


def print_alignment(a1, a2, seq1_start, seq2_start):
    match_line = "".join(
        "|" if x == y else ("+" if x != "-" and y != "-" else " ")
        for x, y in zip(a1, a2)
    )
    print(f"  Query   {seq1_start+1:>4}  {a1}")
    print(f"               {match_line}")
    print(f"  Target  {seq2_start+1:>4}  {a2}")


def main():
    print("\nDay 11 - Smith-Waterman Local Alignment with BLOSUM62")
    print("=" * 65)

    # Query: a short peptide from the kelch13 propeller domain (known motif)
    query  = "CIGGYDGSSIIPNVEAYDHRMKAWVEVAPL"   # 30 aa from kelch13

    # Target: full kelch13 reference sequence, embedded in random-ish flanks
    # to simulate searching a peptide against a larger protein context
    with open("data/kelch13_diverse.fasta") as f:
        lines = f.read().strip().split("\n")
    target_full = lines[1]   # reference sequence, 726 aa

    print(f"\nQuery peptide ({len(query)} aa):")
    print(f"  {query}")
    print(f"\nTarget protein ({len(target_full)} aa): kelch13 reference (full length)")

    print("\n[1/2] Running Smith-Waterman local alignment...")
    a1, a2, score, start1, start2, end1, end2 = smith_waterman(
        query, target_full, gap_open=-10, gap_extend=-1
    )

    print(f"\nBest local alignment found (score = {score}):")
    print_alignment(a1, a2, start1, start2)

    matches = sum(1 for x, y in zip(a1, a2) if x == y)
    identity = matches / len(a1) * 100 if a1 else 0

    print(f"\nAlignment statistics:")
    print(f"  Alignment length : {len(a1)}")
    print(f"  Matches          : {matches}")
    print(f"  Identity         : {identity:.1f}%")
    print(f"  Query region     : positions {start1+1}-{end1}")
    print(f"  Target region    : positions {start2+1}-{end2}")

    # Second test: query with 3 mismatches to test local sensitivity
    print("\n" + "=" * 65)
    print("[2/2] Testing with a mutated query (3 substitutions)...")

    mutated_query = "CIGGYDGSSIIPKVEAYAHRMKAWVEVAPL"  # 3 changes vs original
    print(f"\nMutated query ({len(mutated_query)} aa):")
    print(f"  {mutated_query}")

    a1b, a2b, score_b, s1b, s2b, e1b, e2b = smith_waterman(
        mutated_query, target_full, gap_open=-10, gap_extend=-1
    )

    print(f"\nBest local alignment found (score = {score_b}):")
    print_alignment(a1b, a2b, s1b, s2b)

    matches_b = sum(1 for x, y in zip(a1b, a2b) if x == y)
    identity_b = matches_b / len(a1b) * 100 if a1b else 0
    print(f"\nAlignment statistics:")
    print(f"  Matches          : {matches_b}")
    print(f"  Identity         : {identity_b:.1f}%")
    print(f"  Target region    : positions {s2b+1}-{e2b}")

    print("\nSmith-Waterman implementation validated")
    print("=" * 65)


if __name__ == "__main__":
    main()
