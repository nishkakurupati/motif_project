import motif_finder_utils_noise
import time
import parameters


def find_motifs():
    """Top level function for finding motifs"""
    run_time_stats = []
    for k in range (10):
        for i in range (len(parameters.params_list)):
            ml = parameters.params_list[i][0]
            sc = parameters.params_list[i][1]
            sl = parameters.params_list[i][2]


            start_time = time.perf_counter()
            file_name = str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(k)
            sl,sequences = motif_finder_utils_noise.get_sequences_from_file(file_name)

            motif,positions_of_motif = motif_finder_utils_noise.find_motif(sequences, sl, ml)

            motif_finder_utils_noise.create_output_files(motif, ml,len(sequences),sl, k, positions_of_motif)
            end_time = time.perf_counter()

            # Calculate the duration and print
            runtime = end_time - start_time
            run_time_stats.append([ml, sl, sc, k, runtime])
    return run_time_stats


find_motifs()




