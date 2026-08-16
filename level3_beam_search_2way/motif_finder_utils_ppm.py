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

            
def score_of_candidate(candidate):
    """Sort key: a candidate is [score, motifs, sites], rank by score"""
    return candidate[0]


def keep_top_n(candidates, beam_width):
    """Sort candidates by score (highest first) and keep the top beam_width.
    Python's sort is stable, so equal scores keep first-seen order (first wins)"""
    candidates.sort(key=score_of_candidate, reverse=True)
    return candidates[:beam_width]


def run_one_trial(sequences, sl, ml, beam_width):
    """One randomized beam-search trial: shuffle the sequence order, seed from
    all pairs of the first two sequences keeping the top beam_width partial
    motifs, then for each remaining sequence extend every survivor by every
    candidate l-mer and keep the top beam_width again"""

    # Random order of sequences (tracked by original index)
    order = []
    for i in range (len(sequences)):
        order.append(i)
    random.shuffle(order)

    # Seed: all pairs of seq order[0] x order[1], keep top beam_width.
    # Each candidate is [score, motifs, sites] with sites indexed by ORIGINAL index.
    seq1_motifs = get_possible_motifs(sl, ml, sequences[order[0]])
    seq2_motifs = get_possible_motifs(sl, ml, sequences[order[1]])

    candidates = []
    for i in range (len(seq1_motifs)):
        for j in range (len(seq2_motifs)):
            motifs = [seq1_motifs[i], seq2_motifs[j]]
            score = get_score(build_ppm_from_brute_force(motifs))

            sites = []
            for s in range (len(sequences)):
                sites.append(0)
            sites[order[0]] = i
            sites[order[1]] = j

            candidates.append([score, motifs, sites])

    beam = keep_top_n(candidates, beam_width)

    # Extend: for each remaining sequence, grow every beam entry by every l-mer
    for t in range(2, len(sequences)):
        seq_index = order[t]
        possible_motifs = get_possible_motifs(sl, ml, sequences[seq_index])

        new_candidates = []
        for b in range (len(beam)):
            entry_motifs = beam[b][1]
            entry_sites = beam[b][2]
            for k in range (len(possible_motifs)):
                test_motifs = entry_motifs + [possible_motifs[k]]
                new_score = get_score(build_ppm_from_brute_force(test_motifs))

                new_sites = []
                for s in range (len(entry_sites)):
                    new_sites.append(entry_sites[s])
                new_sites[seq_index] = k

                new_candidates.append([new_score, test_motifs, new_sites])

        beam = keep_top_n(new_candidates, beam_width)

    best = beam[0]
    best_score = best[0]
    best_motifs = best[1]
    best_sites = best[2]
    best_ppm = build_ppm_from_brute_force(best_motifs)

    return best_motifs, best_ppm, best_sites, best_score


def greedy_motif_search(sequences, sl, ml, num_trials, beam_width):
    """Run multiple randomized order trials and keep the one with the best PPM score"""
    best_trial_score = -1
    best_trial_motifs = []
    best_trial_ppm = None
    best_trial_sites = []

    for trial in range (num_trials):
        motifs, ppm, sites, score = run_one_trial(sequences, sl, ml, beam_width)

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



