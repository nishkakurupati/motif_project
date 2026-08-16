"""Fast-iteration harness: run the engine on all 7 param buckets over a small
k-range, then report sites-matched (the key metric) vs the best incumbent,
plus KL divergence and runtime.

Usage:
    python run_subset.py                 # all buckets, k=0..2, select=profile
    python run_subset.py 5               # k=0..4
    python run_subset.py 3 logodds       # k=0..2, log-odds selection
"""
import sys
import time
import random

import motif_finder_utils_ppm
import motif_benchmark_ppm
import parameters

# Best sites-matched among the old variants (randomized / topN / 6way) at p=0.8,
# keyed by (ml, sc, sl). This is the bar to beat.
INCUMBENT = {
    (6, 10, 500): 0.05,
    (7, 10, 500): 0.19,
    (8, 10, 500): 0.19,
    (8, 5, 500): 0.08,
    (8, 20, 500): 0.31,
    (8, 10, 1000): 0.06,
    (8, 10, 2000): 0.00,
}


def run(num_k, select):
    print(f"knobs: restarts=adaptive sweeps={parameters.gibbs_sweeps} "
          f"alpha={parameters.pseudocount} em={parameters.em_max_iters} select={select} k=0..{num_k-1}")
    print(f"{'config':>12} {'sites':>7} {'incumb':>7} {'delta':>7} {'KL':>8} {'sec/ds':>7}")

    wins = 0
    for cfg in parameters.params_list:
        ml, sc, sl = cfg
        random.seed(1)
        matched = 0
        total = 0
        kl_sum = 0.0
        t0 = time.perf_counter()
        for k in range(num_k):
            file_name = str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(k)
            _sl, sequences = motif_finder_utils_ppm.get_sequences_from_file(file_name)
            best_ppm, sites = motif_finder_utils_ppm.find_motif_sites(
                sequences, _sl, ml,
                parameters.restarts_for(_sl), parameters.gibbs_sweeps,
                parameters.pseudocount, parameters.em_max_iters,
                0.8, select,
            )
            motif_finder_utils_ppm.create_output_files(best_ppm, ml, len(sequences), _sl, k, sites)
            _, count = motif_benchmark_ppm.compare_sites(ml, sc, sl, k)
            kl_sum += motif_benchmark_ppm.compare_ppm(ml, sc, sl, k)
            matched += count
            total += sc
        elapsed = time.perf_counter() - t0
        sites_frac = matched / total
        inc = INCUMBENT.get((ml, sc, sl), 0.0)
        delta = sites_frac - inc
        if delta > 0:
            wins += 1
        flag = "  BEAT" if delta > 0 else ("  tie" if delta == 0 else "")
        print(f"{ml}_{sc}_{sl:<4} {sites_frac:>7.3f} {inc:>7.3f} {delta:>+7.3f} "
              f"{kl_sum/num_k:>8.3f} {elapsed/num_k:>7.2f}{flag}")
    print(f"buckets beaten: {wins}/{len(parameters.params_list)}")


if __name__ == "__main__":
    num_k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    select = sys.argv[2] if len(sys.argv) > 2 else 'profile'
    run(num_k, select)
