# scripts/day11_blosum_comparison.py
# Author: Boris Djagou
# Date: July 15, 2026
# Compare scoring behavior across different BLOSUM matrices

from Bio.Align import substitution_matrices

def main():
    print("\nDay 11 - BLOSUM Matrix Comparison")
    matrices = {
        "BLOSUM45": substitution_matrices.load("BLOSUM45"),
        "BLOSUM62": substitution_matrices.load("BLOSUM62"),
        "BLOSUM80": substitution_matrices.load("BLOSUM80"),
    }

    # Test pairs: conservative vs non-conservative substitutions
    pairs = [
        ("L", "I", "Leu -> Ile (conservative, both hydrophobic)"),
        ("K", "R", "Lys -> Arg (conservative, both basic)"),
        ("G", "W", "Gly -> Trp (non-conservative, size mismatch)"),
        ("D", "E", "Asp -> Glu (conservative, both acidic)"),
        ("C", "Y", "Cys -> Tyr (kelch13 C580Y resistance mutation)"),
        ("A", "A", "Ala -> Ala (identity)"),
    ]

    print(f"\n{'Substitution':<45} {'BLOSUM45':>10} {'BLOSUM62':>10} {'BLOSUM80':>10}")
    for aa1, aa2, description in pairs:
        scores = []
        for name, matrix in matrices.items():
            try:
                s = matrix[aa1, aa2]
            except KeyError:
                s = matrix[aa2, aa1]
            scores.append(s)
        print(f"{description:<45} {scores[0]:>10} {scores[1]:>10} {scores[2]:>10}")

    print("\nInterpretation:")
    print("  Higher score = more likely substitution in evolution")
    print("  Negative score = unlikely / destabilizing substitution")
    print("  BLOSUM80 gives harsher penalties (for close homologs)")
    print("  BLOSUM45 is more permissive (for distant homologs)")

    print("\nDone")


if __name__ == "__main__":
    main()
