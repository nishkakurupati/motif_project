import check_accuracy_ppm
import motif_finder_ppm
import evaluate_ppm_utils

run_time_stats = motif_finder_ppm.find_motifs()

motif_match_stats, sites_match_stats = check_accuracy_ppm.check_acc()

print("--------Runtime Stats-----------")
evaluate_ppm_utils.evaluate(run_time_stats)
print("--------Accuracy Stats-----------")
evaluate_ppm_utils.evaluate(motif_match_stats)
print("--------Motif Sites Stats-----------")
evaluate_ppm_utils.evaluate(sites_match_stats, True)
