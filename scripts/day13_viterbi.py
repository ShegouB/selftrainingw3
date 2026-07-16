# scripts/day13_viterbi.py
# Author: Boris Djagou
# Date: July 17, 2026
# Simplified Viterbi algorithm to understand HMM decoding

import numpy as np


def viterbi(observations, states, start_prob, trans_prob, emit_prob):
    """
    Standard Viterbi algorithm implementation.
    observations : list of observed symbols
    states       : list of hidden state names
    start_prob   : dict of initial state probabilities
    trans_prob   : dict of dict, transition probabilities between states
    emit_prob    : dict of dict, emission probabilities per state
    Returns: best path (list of states), probability of that path
    """
    n_obs = len(observations)
    n_states = len(states)

    V = [{}]
    path = {}

    for state in states:
        V[0][state] = start_prob[state] * emit_prob[state][observations[0]]
        path[state] = [state]

    for t in range(1, n_obs):
        V.append({})
        new_path = {}

        for curr_state in states:
            best_prob, best_prev = max(
                (V[t - 1][prev_state] * trans_prob[prev_state][curr_state] *
                 emit_prob[curr_state][observations[t]], prev_state)
                for prev_state in states
            )
            V[t][curr_state] = best_prob
            new_path[curr_state] = path[best_prev] + [curr_state]

        path = new_path

    best_final_prob, best_final_state = max(
        (V[n_obs - 1][state], state) for state in states
    )

    return path[best_final_state], best_final_prob


def demo_cpg_island_detector():
    """
    Classic bioinformatics HMM example: detect CpG islands in DNA.
    Two hidden states: 'CpG' (high GC region) vs 'normal' (background).
    This mirrors detecting AT-rich vs GC-rich regions, directly
    relevant to Project 1's genomic island detection.
    """
    print("\nDemo: CpG Island Detection via Viterbi (2-state HMM)")
    print("-" * 60)

    states = ["CpG", "normal"]
    observations = list("GCGCGATATATGCGCATATATATGCGCGC")

    start_prob = {"CpG": 0.5, "normal": 0.5}

    trans_prob = {
        "CpG":    {"CpG": 0.85, "normal": 0.15},
        "normal": {"CpG": 0.10, "normal": 0.90},
    }

    emit_prob = {
        "CpG":    {"A": 0.10, "T": 0.10, "G": 0.40, "C": 0.40},
        "normal": {"A": 0.30, "T": 0.30, "G": 0.20, "C": 0.20},
    }

    best_path, best_prob = viterbi(observations, states, start_prob,
                                    trans_prob, emit_prob)

    print(f"\nSequence:  {''.join(observations)}")
    print(f"Best path: {''.join('C' if s == 'CpG' else 'n' for s in best_path)}")
    print(f"           (C = predicted CpG/GC-rich island, n = normal)")
    print(f"\nLog probability of best path: {np.log(best_prob):.4f}")

    cpg_regions = []
    in_cpg = False
    start = 0
    for i, s in enumerate(best_path):
        if s == "CpG" and not in_cpg:
            start = i
            in_cpg = True
        elif s != "CpG" and in_cpg:
            cpg_regions.append((start, i))
            in_cpg = False
    if in_cpg:
        cpg_regions.append((start, len(best_path)))

    print(f"\nPredicted CpG/GC-rich islands:")
    for start, end in cpg_regions:
        segment = "".join(observations[start:end])
        print(f"  Position {start}-{end}: {segment}")


def main():
    print("\nDay 13 - Hidden Markov Models and the Viterbi Algorithm")
    print("=" * 65)

    demo_cpg_island_detector()

    print("\n" + "=" * 65)
    print("Connection to Week 1 Project 1:")
    print("  Your GC-content scanner used a SLIDING WINDOW approach")
    print("  to find AT-rich regions in P. falciparum.")
    print("  An HMM-based approach (like this Viterbi demo) would")
    print("  instead PROBABILISTICALLY SEGMENT the genome into GC-rich")
    print("  and AT-rich states, which is exactly how tools like")
    print("  GeneMark and Augustus detect gene boundaries in genomes.")

    print("\nDone")
    print("=" * 65)


if __name__ == "__main__":
    main()
