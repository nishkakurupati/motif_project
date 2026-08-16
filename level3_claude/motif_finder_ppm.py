import motif_finder_utils_ppm
import time
import random
import parameters


def find_motifs():
    """Top level function for finding motifs (Gibbs + hard-EM engine)."""
    # Fixed seed so the randomized search is reproducible run-to-run.
    random.seed(1)

    run_time_stats = []
    for k in range(10):
        for i in range(len(parameters.params_list)):
            ml = parameters.params_list[i][0]
            sc = parameters.params_list[i][1]
            sl = parameters.params_list[i][2]
            print(f"Running {i} of {len(parameters.params_list)}, k={k}  (ml={ml}, sc={sc}, sl={sl})")

            start_time = time.perf_counter()
            file_name = str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(k)
            sl, sequences = motif_finder_utils_ppm.get_sequences_from_file(file_name)

            best_ppm, sites = motif_finder_utils_ppm.find_motif_sites(
                sequences, sl, ml,
                parameters.restarts_for(sl),
                parameters.gibbs_sweeps,
                parameters.pseudocount,
                parameters.em_max_iters,
                0.8, 'logodds',
            )

            motif_finder_utils_ppm.create_output_files(best_ppm, ml, len(sequences), sl, k, sites)
            end_time = time.perf_counter()

            runtime = end_time - start_time
            run_time_stats.append([ml, sl, sc, k, runtime])
    return run_time_stats


# Only auto-run when executed directly (python motif_finder_ppm.py).
# evaluate_ppm.py imports this module and calls find_motifs() itself.
if __name__ == "__main__":
    find_motifs()
