import check_accuracy
import motif_finder

def evaluate(stat_stats, use_sites = False):
    """Evaluate performance of motif finders"""
    evaluate_ml(stat_stats, use_sites)
    evaluate_sl(stat_stats, use_sites)
    evaluate_sc(stat_stats, use_sites)


def evaluate_ml(stats_object, use_sites = False):
    """Evaluating the performance of changing motif length"""
    ml_6_stats = []
    ml_7_stats = []
    ml_8_stats = []

    for i in range (len(stats_object)):
        ml,sl,sc, index, stat = stats_object[i] 
        if sl == 500 and sc == 10:
            if ml == 6:
                ml_6_stats.append(stat)
            elif ml == 7:
                ml_7_stats.append(stat)
            else:
                ml_8_stats.append(stat)           

    if use_sites:
        total = len(ml_6_stats) * 10
    else:
        total = len(ml_6_stats)  # same for ml_7 and 8

    ml_6_stat = sum(ml_6_stats)/total
    ml_7_stat = sum(ml_7_stats)/total
    ml_8_stat = sum(ml_8_stats)/total

    print("ml 6 average stats ", ml_6_stat)
    print("ml 7 average stats ", ml_7_stat)
    print("ml 8 average stats ", ml_8_stat)



def evaluate_sl(stats_object, use_sites = False):
    """Evaluating the performance of changing sequence length"""
    sl_500_stats = []
    sl_1000_stats = []
    sl_2000_stats = []


    for i in range (len(stats_object)):
        ml,sl,sc, index, stat = stats_object[i]
        if ml == 8 and sc == 10:
            if sl == 500:
                sl_500_stats.append(stat)
            elif sl == 1000:
                sl_1000_stats.append(stat)
            else:
                sl_2000_stats.append(stat)

        if use_sites:
            total = len(sl_500_stats) * 10
        else:
            total = len(sl_500_stats)  # same for sl_1000 and 2000
            


    sl_500_stat = sum(sl_500_stats)/total 
    sl_1000_stat = sum(sl_1000_stats)/total
    sl_2000_stat = sum(sl_2000_stats)/total

    print("sl 500 average stats ", sl_500_stat)
    print("sl 1000 average stats ", sl_1000_stat)
    print("sl 2000 average stats ", sl_2000_stat)


def evaluate_sc(stats_object, use_sites = False):
    """Evaluating the performance of changing sequence count"""
    sc_5_stats = []
    sc_10_stats = []
    sc_20_stats = []

    total_sites = 0
    for i in range (len(stats_object)):
        ml,sl,sc, index, stat = stats_object[i]
        if ml == 8 and sl == 500:
            if sc == 5:
                sc_5_stats.append(stat)
                total_sites += sc
            elif sc == 10:
                sc_10_stats.append(stat)
            else:
                sc_20_stats.append(stat)
            total_sites += sc

        if use_sites:
            sc_5_total = len(sc_5_stats) * 5
            sc_10_total  = len(sc_10_stats) * 10
            sc_20_total = len(sc_20_stats) * 20
        else:
            sc_5_total = len(sc_5_stats) 
            sc_10_total  = len(sc_10_stats)
            sc_20_total = len(sc_20_stats)


    sc_5_stat = sum(sc_5_stats)/sc_5_total 
    sc_10_stat = sum(sc_10_stats)/sc_10_total 
    sc_20_stat = sum(sc_20_stats)/sc_20_total 

    print("sc 5 average stats ", sc_5_stat)
    print("sc 10 average stats ", sc_10_stat)
    print("sc 20 average stats ", sc_20_stat)

