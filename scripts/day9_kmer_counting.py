# scripts/day9_kmer_counting.py
# Author: Boris Djagou
# Date: July 13-14, 2026 (catch-up)
# Exercise: Count all 4-mers in a 1000-nt sequence

from collections import defaultdict
import random

def generate_random_sequence(length, seed=42):
    """Generate a random DNA sequence for testing."""
    random.seed(seed)
    return "".join(random.choice("ACGT") for _ in range(length))


def count_kmers(sequence, k):
    """Count all k-mers in a sequence using a sliding window."""
    kmer_counts = defaultdict(int)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] += 1
    return dict(kmer_counts)


def main():
    print("\n Day 9 catch-up — k-mer Counting")

    # Use your real kelch13 reference sequence if available, else random
    try:
        with open("data/kelch13_diverse.fasta") as f:
            lines = f.read().strip().split("\n")
        # Get first sequence (reference)
        seq = lines[1]
        print(f"Using real kelch13 reference sequence ({len(seq)} aa)")
    except FileNotFoundError:
        seq = generate_random_sequence(1000)
        print(f"Using random 1000-nt test sequence")

    k = 4
    kmers = count_kmers(seq, k)

    print(f"\nSequence length: {len(seq)}")
    print(f"Total possible {k}-mers: {len(seq) - k + 1}")
    print(f"Unique {k}-mers found: {len(kmers)}")

    # Show top 10 most frequent k-mers
    sorted_kmers = sorted(kmers.items(), key=lambda x: x[1], reverse=True)
    print(f"\nTop 10 most frequent {k}-mers:")
    print(f"  {'k-mer':<8} {'Count':>6}")
    for kmer, count in sorted_kmers[:10]:
        print(f"  {kmer:<8} {count:>6}")

    print(f"\n k-mer counting complete")


if __name__ == "__main__":
    main()
