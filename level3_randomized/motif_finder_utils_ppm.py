import random
import math

def get_sequences_from_file(file_name):
    sequences = []

    filename = 'data_set/file_' + file_name + "sequences.fa"
    file_obj = open(filename, 'r')
    for line in file_obj:
        sequences.append(line.rstrip('\n'))
        sl = len(line.rstrip('\n'))
    file_obj.close()

    return sl, sequences

def get_possible_motifs(sl,ml, sequence):
    possible_motifs = []

    for i in range (sl - ml + 1):
        pos_motif = sequence[i:i + ml]

        possible_motifs.append(pos_motif)

    return possible_motifs

def build_ppm_from_brute_force(motifs):
    ppm = []
    for i in range (len(motifs[0])):
        dict_counter = dict()
        dict_counter['A'] = 0
        dict_counter['T'] = 0
        dict_counter['G'] = 0
        dict_counter['C'] = 0
        for motif in motifs:
            dict_counter[motif[i]] += 1

        for key in dict_counter:
            dict_counter[key] = dict_counter[key]/len(motifs)

        ppm.append(dict_counter)
    return ppm
        

def get_score(ppm):
    total_score = 0
    for column in ppm:
        score = 0
        for probability in column.values():
            if probability > 0:
                score += probability * math.log2(4 * probability)
        total_score += score

    return total_score


def compare_ppms(ppms):
    best_ppm = []
    best_score = -1
    for ppm in ppms:
        score = get_score(ppm)
        if score > best_score:
            best_score = score
            best_ppm = ppm

    return best_ppm, best_score

            
def run_one_trial(sequences, sl, ml):
    """One randomized greedy trial: shuffle the sequence order, all-pairs seed
    the first two, then greedily add the rest in that random order"""

    # Random order of sequences (tracked by original index)
    order = []
    for i in range (len(sequences)):
        order.append(i)
    random.shuffle(order)

    # sites is indexed by the ORIGINAL sequence index
    sites = []
    for i in range (len(sequences)):
        sites.append(0)

    # All-pairs seed on the first two sequences of the random order
    seq1_motifs = get_possible_motifs(sl, ml, sequences[order[0]])
    seq2_motifs = get_possible_motifs(sl, ml, sequences[order[1]])
    best_score = -1
    best_motifs = []
    best_i = 0
    best_j = 0

    for i in range (len(seq1_motifs)):
        for j in range (len(seq2_motifs)):
            current_motifs = [seq1_motifs[i], seq2_motifs[j]]

            current_ppm = build_ppm_from_brute_force(current_motifs)
            score = get_score(current_ppm)

            if score > best_score:
                best_score = score
                best_motifs = current_motifs
                best_i, best_j = i,j

    sites[order[0]] = best_i
    sites[order[1]] = best_j

    # Greedily add the remaining sequences in the random order
    for t in range(2, len(sequences)):
        seq_index = order[t]
        possible_motifs = get_possible_motifs(sl, ml, sequences[seq_index])
        best_local_score = -1
        best_local_motif = None
        best_k = 0
        for k in range (len(possible_motifs)):
            motif = possible_motifs[k]
            test_motifs = best_motifs + [motif]

            ppm = build_ppm_from_brute_force(test_motifs)
            score = get_score(ppm)

            if score > best_local_score:
                best_local_score = score
                best_local_motif = motif
                best_k = k
        sites[seq_index] = best_k
        best_motifs.append(best_local_motif)

    final_ppm = build_ppm_from_brute_force(best_motifs)
    final_score = get_score(final_ppm)

    return best_motifs, final_ppm, sites, final_score


def greedy_motif_search(sequences, sl, ml, num_trials):
    """Run multiple randomized trials and keep the one with the best PPM score"""
    best_trial_score = -1
    best_trial_motifs = []
    best_trial_ppm = None
    best_trial_sites = []

    for trial in range (num_trials):
        motifs, ppm, sites, score = run_one_trial(sequences, sl, ml)

        if score > best_trial_score:
            best_trial_score = score
            best_trial_motifs = motifs
            best_trial_ppm = ppm
            best_trial_sites = sites

    return best_trial_motifs, best_trial_ppm, best_trial_sites

def create_output_files(ppm, ml, sc, sl, index, sites):
    start_name = "data_set/file_" + str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(index)
    filename = start_name + "predicted_motif.txt"

    file_obj = open(filename, 'w')

    for n in ['A', 'T', 'C', 'G']:
        row = ""
        for col in range(len(ppm)):
            row += str(ppm[col][n]) + " "

        file_obj.write(n + " " + row.strip() + "\n")
    file_obj.close()

    filename2 = start_name + "predicted_sites.txt"

    file_obj = open(filename2, 'w')
    for i in range (sc):
        file_obj.write(str(sites[i]) + "\n")
    file_obj.close()



