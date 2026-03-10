import check_accuracy_with_noise
import motif_finder_noise
import evaluate_utils_with_noise

run_time_stats = motif_finder_noise.find_motifs()

motif_match_stats, sites_match_stats = check_accuracy_with_noise.check_acc()

print("--------Runtime Stats-----------")
evaluate_utils_with_noise.evaluate(run_time_stats)
print("--------Accuracy Stats-----------")
evaluate_utils_with_noise.evaluate(motif_match_stats)
print("--------Motif Sites Stats-----------")
evaluate_utils_with_noise.evaluate(sites_match_stats,True)
