def compare_motif(ml,sc,sl, index):
    """Comparing if the predicted motif matched the actual motif"""
    is_same = False

    start_name = "data_set/file_" + str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(index)
    filename = start_name + "motif.txt"
    file_obj = open(filename, 'r')
    for line in file_obj:
	    motif = line
    file_obj.close()

    filename_pred = start_name + "predictedmotif.txt"
    file_obj = open(filename_pred, 'r')
    for line in file_obj:
	    motif_pred = line
    file_obj.close()

    if motif == motif_pred:
        is_same = True

    return is_same


def compare_sites(ml,sc,sl, index):
    """Comparing if the predicted sites matched the actual sites"""
    is_same = False

    sites = []
    pred_sites = []
    start_name = "data_set/file_" + str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(index)

    filename = start_name + "sites.txt"
    file_obj = open(filename, 'r')
    for line in file_obj:
	    sites.append(line)
    file_obj.close()

    filename_pred = start_name + "predictedsites.txt"
    file_obj = open(filename_pred, 'r')
    for line in file_obj:
	    pred_sites.append(line)
    file_obj.close()

    count = 0
    for i in range (sc):
        if sites[i] == pred_sites[i]:
            count = count + 1
    
    if count == sc:
        is_same = True

    return is_same, count
