import check_accuracy
import motif_finder
import evaluate_utils

run_time_stats = motif_finder.find_motifs()

motif_match_stats, sites_match_stats = check_accuracy.check_acc()

print("--------Runtime Stats-----------")
evaluate_utils.evaluate(run_time_stats)
print("--------Accuracy Stats-----------")
evaluate_utils.evaluate(motif_match_stats)
print("--------Motif Sites Stats-----------")
evaluate_utils.evaluate(sites_match_stats,True)

