# scripts/day12_blast_from_scratch.py
# Author: Boris Djagou
# Date: July 16, 2026
# Simplified implementation of BLAST word-matching and extension heuristics

from Bio.Align import substitution_matrices
import math

BLOSUM62 = substitution_matrices.load("BLOSUM62")


def get_score(a, b, matrix=BLOSUM62):
    try:
        return matrix[a, b]
    except KeyError:
        return matrix[b, a]


def build_word_index(sequence, word_size=3):
    """Build an index of all words (k-mers) and their positions in a sequence."""
    index = {}
    for i in range(len(sequence) - word_size + 1):
        word = sequence[i:i + word_size]
        index.setdefault(word, []).append(i)
    return index


def find_seed_matches(query, target, word_size=3):
    """Step 1: find exact word matches (seeds) between query and target."""
    query_words = build_word_index(query, word_size)
    target_words = build_word_index(target, word_size)

    seeds = []
    for word, query_positions in query_words.items():
        if word in target_words:
            for qpos in query_positions:
                for tpos in target_words[word]:
                    seeds.append((qpos, tpos, word))
    return seeds


def extend_seed(query, target, qpos, tpos, word_size=3, drop_threshold=5):
    """Step 2: extend a seed match in both directions using BLOSUM62 scoring.
    Stop extension when score drops drop_threshold below the best score seen.
    """
    score = sum(get_score(query[qpos + i], target[tpos + i]) for i in range(word_size))
    best_score = score

    left_ext = 0
    right_ext = 0

    # Extend right
    running_score = score
    i = word_size
    while qpos + i < len(query) and tpos + i < len(target):
        running_score += get_score(query[qpos + i], target[tpos + i])
        if running_score > best_score:
            best_score = running_score
            right_ext = i - word_size + 1
        if best_score - running_score > drop_threshold:
            break
        i += 1

    # Extend left
    running_score = best_score
    i = 1
    while qpos - i >= 0 and tpos - i >= 0:
        running_score += get_score(query[qpos - i], target[tpos - i])
        if running_score > best_score:
            best_score = running_score
            left_ext = i
        if best_score - running_score > drop_threshold:
            break
        i += 1

    start_q = qpos - left_ext
    start_t = tpos - left_ext
    length = left_ext + word_size + right_ext

    return {
        "score": best_score,
        "query_start": start_q,
        "target_start": start_t,
        "length": length,
        "query_seq": query[start_q:start_q + length],
        "target_seq": target[start_t:start_t + length],
    }


def compute_evalue(score, query_len, db_len, lambda_param=0.267, k_param=0.041):
    """Simplified Karlin-Altschul E-value estimate.
    E = K * m * n * e^(-lambda * S)
    """
    E = k_param * query_len * db_len * math.exp(-lambda_param * score)
    return E


# In simplified_blast(), replace the seen_regions check with an HSP-level dedup:

def simplified_blast(query, target, word_size=3, drop_threshold=5, min_score=20):
    seeds = find_seed_matches(query, target, word_size)

    hsps = []
    seen_hsps = set()

    for qpos, tpos, word in seeds:
        hsp = extend_seed(query, target, qpos, tpos, word_size, drop_threshold)
        hsp_key = (hsp["query_start"], hsp["target_start"], hsp["length"])

        if hsp["score"] >= min_score and hsp_key not in seen_hsps:
            hsps.append(hsp)
            seen_hsps.add(hsp_key)

    hsps.sort(key=lambda x: x["score"], reverse=True)
    return hsps, len(seeds)


def main():
    print("\nDay 12 - Simplified BLAST Implementation (word matching + extension)")

    mutated_query = "CIGGYDGSSIIPKVEAYAHRMKAWVEVAPL"

    with open("data/kelch13_diverse.fasta") as f:
        content = f.read().strip()

    entries = content.split(">")[1:]
    target = None
    for entry in entries:
        lines = entry.split("\n")
        if "REFERENCE" in lines[0]:
            target = "".join(lines[1:])
            break

    print(f"\nQuery ({len(mutated_query)} aa): {mutated_query}")
    print(f"Target ({len(target)} aa): kelch13 reference")

    print("\n[1/3] Building word index and finding seed matches (word size = 3)...")
    seeds = find_seed_matches(mutated_query, target, word_size=3)
    print(f"  Total seed matches found: {len(seeds)}")

    print("\n[2/3] Extending seeds and filtering by minimum score...")
    hsps, n_seeds = simplified_blast(mutated_query, target, word_size=3,
                                      drop_threshold=5, min_score=20)
    print(f"  High-scoring segment pairs (HSPs) after extension: {len(hsps)}")

    print("\n[3/3] Top hits with E-value estimation:")
    db_len = len(target)

    for i, hsp in enumerate(hsps[:5]):
        evalue = compute_evalue(hsp["score"], len(mutated_query), db_len)
        print(f"\n  Hit {i+1}: score={hsp['score']}, length={hsp['length']}, "
              f"E-value={evalue:.2e}")
        print(f"    Query  {hsp['query_start']+1:>4}  {hsp['query_seq']}")
        print(f"    Target {hsp['target_start']+1:>4}  {hsp['target_seq']}")

    print("\nSimplified BLAST implementation complete")


if __name__ == "__main__":
    main()
