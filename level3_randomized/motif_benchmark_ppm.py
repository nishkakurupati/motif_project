import math


def compare_ppm(ml,sc,sl, index):
    """Comparing predicted PPM against actual PPM using KL divergence"""
    start_name = "data_set/file_" + str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(index)

    filename = start_name + "motif.txt"
    file_obj = open(filename, 'r')
    ppm_rows = dict()
    for line in file_obj:
        parts = line.rstrip('\n').split()
        nucleotide = parts[0]
        values = []
        for j in range(1, len(parts)):
            values.append(float(parts[j]))
        ppm_rows[nucleotide] = values
    file_obj.close()

    ml_len = len(ppm_rows['A'])
    actual_ppm = []
    for i in range(ml_len):
        column = dict()
        for n in ['A', 'T', 'C', 'G']:
            column[n] = ppm_rows[n][i]
        actual_ppm.append(column)

    filename_pred = start_name + "predicted_motif.txt"
    file_obj = open(filename_pred, 'r')
    pred_rows = dict()
    for line in file_obj:
        parts = line.rstrip('\n').split()
        nucleotide = parts[0]
        values = []
        for j in range(1, len(parts)):
            values.append(float(parts[j]))
        pred_rows[nucleotide] = values
    file_obj.close()

    predicted_ppm = []
    for i in range(ml_len):
        column = dict()
        for n in ['A', 'T', 'C', 'G']:
            column[n] = pred_rows[n][i]
        predicted_ppm.append(column)

    total_score = 0
    for i in range(len(predicted_ppm)):
        predicted_column = predicted_ppm[i]
        actual_column = actual_ppm[i]
        nucleotides = ['A', 'C', 'G', 'T']
        for nucleotide in nucleotides:
            predicted_value = predicted_column[nucleotide]
            actual_value = actual_column[nucleotide]
            if predicted_value == 0:
                continue
            if actual_value == 0:
                actual_value = 0.0001
            value = predicted_value * math.log(predicted_value / actual_value, 2)
            total_score = total_score + value

    return total_score


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

    filename_pred = start_name + "predicted_sites.txt"
    file_obj = open(filename_pred, 'r')
    for line in file_obj:
        pred_sites.append(line)
    file_obj.close()

    count = 0
    for i in range(sc):
        if sites[i] == pred_sites[i]:
            count = count + 1

    if count == sc:
        is_same = True

    return is_same, count
