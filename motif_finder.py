import motif_finder_utils
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
            ml = motif_finder_utils.get_ml_from_file(file_name)
            sl,sequences = motif_finder_utils.get_sequences_from_file(file_name)
        
            possible_motifs = motif_finder_utils.create_possible_motifs(sl,ml,sequences)
            possible_motifs = motif_finder_utils.find_motif(possible_motifs, sequences, sl, ml)
            positions_of_motif = motif_finder_utils.get_sites(sl,ml,possible_motifs, sequences)

            motif_finder_utils.create_output_files(possible_motifs, ml,len(sequences),sl, k, positions_of_motif)
            end_time = time.perf_counter()

            # Calculate the duration and print
            runtime = end_time - start_time
            run_time_stats.append([ml, sl, sc, k, runtime])
    return run_time_stats






