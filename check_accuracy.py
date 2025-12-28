import motif_benchmark
import parameters


def check_acc():
    motif_match_stats = []
    sites_match_stats = []

    for k in range (10):
        for i in range (len(parameters.params_list)):
            ml = parameters.params_list[i][0]
            sc = parameters.params_list[i][1]
            sl = parameters.params_list[i][2]

            motif_match = motif_benchmark.compare_motif(ml,sc,sl,k)
            
            full_site_match, site_matches = motif_benchmark.compare_sites(ml,sc,sl,k)



            motif_match_stats.append([ml, sl, sc, k, motif_match])
            sites_match_stats.append([ml, sl, sc, k, site_matches])

    return motif_match_stats, sites_match_stats