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

            
def greedy_motif_search(sequences, sl, ml):
    best_motifs = [] 
    best_ppm = None 
    best_score = -1 
    sites = []

    seq1_motifs = get_possible_motifs(sl, ml, sequences[0]) 
    seq2_motifs = get_possible_motifs(sl, ml, sequences[1]) 
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

    sites.append(best_i)
    sites.append(best_j)

    for i in range(2, len(sequences)):
        possible_motifs = get_possible_motifs(sl, ml, sequences[i])
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
        sites.append(best_k)
        best_motifs.append(best_local_motif)

    best_ppm = build_ppm_from_brute_force(best_motifs)
    
    return best_motifs, best_ppm, sites

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



