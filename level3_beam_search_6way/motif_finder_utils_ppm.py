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


def seed_pair_topN(sequences, sl, ml, a, b, beam_width):
    """Exhaustive all-pairs seed on two sequences (original indices a, b).
    Returns the top beam_width candidates, each [score, motifs, sites] with
    sites a full-length array (indexed by ORIGINAL index) with only a and b set"""
    seqa_motifs = get_possible_motifs(sl, ml, sequences[a])
    seqb_motifs = get_possible_motifs(sl, ml, sequences[b])

    candidates = []
    for i in range (len(seqa_motifs)):
        for j in range (len(seqb_motifs)):
            motifs = [seqa_motifs[i], seqb_motifs[j]]
            score = get_score(build_ppm_from_brute_force(motifs))

            sites = []
            for s in range (len(sequences)):
                sites.append(0)
            sites[a] = i
            sites[b] = j

            candidates.append([score, motifs, sites])

    return keep_top_n(candidates, beam_width)


def single_sequence_fallback(sequences, sl, ml, beam_width, order):
    """Degenerate case: only one sequence to seed from. Beam is the top
    beam_width l-mers of that sequence, each a single-motif candidate"""
    seq_index = order[0]
    possible_motifs = get_possible_motifs(sl, ml, sequences[seq_index])

    candidates = []
    for k in range (len(possible_motifs)):
        motifs = [possible_motifs[k]]
        score = get_score(build_ppm_from_brute_force(motifs))

        sites = []
        for s in range (len(sequences)):
            sites.append(0)
        sites[seq_index] = k

        candidates.append([score, motifs, sites])

    beam = keep_top_n(candidates, beam_width)
    best = beam[0]
    best_ppm = build_ppm_from_brute_force(best[1])
    return best[1], best_ppm, best[2], best[0]


def run_one_trial(sequences, sl, ml, beam_width):
    """One randomized beam-search trial with multi-group seeding: shuffle the
    sequence order, split the first sequences into up to 3 disjoint groups of 2,
    seed each group exhaustively (own top beam_width), join all groups by
    cross-product (keep top beam_width), then extend the remaining sequences
    sequentially like the topN variant"""

    # Random order of sequences (tracked by original index)
    order = []
    for i in range (len(sequences)):
        order.append(i)
    random.shuffle(order)

    t = len(sequences)
    # As many disjoint pairs as we have, capped at 3 groups
    num_groups = min(3, t // 2)

    # Degenerate: 0 or 1 sequence, no pair to seed
    if num_groups == 0:
        return single_sequence_fallback(sequences, sl, ml, beam_width, order)

    # Per-group exhaustive all-pairs seed, each keeps its own top beam_width.
    # Remember each group's two original indices for a safe site merge.
    group_beams = []
    group_pairs = []
    for g in range (num_groups):
        a = order[2 * g]
        b = order[2 * g + 1]
        group_beams.append(seed_pair_topN(sequences, sl, ml, a, b, beam_width))
        group_pairs.append([a, b])

    # Join groups: cross-product, merge motifs and sites, keep top beam_width.
    # Merge copies group g's site values only at its own known indices (a, b),
    # so a real motif at position 0 is never mistaken for an empty slot.
    beam = group_beams[0]
    for g in range (1, num_groups):
        a = group_pairs[g][0]
        b = group_pairs[g][1]
        joined = []
        for x in range (len(beam)):
            x_motifs = beam[x][1]
            x_sites = beam[x][2]
            for y in range (len(group_beams[g])):
                y_motifs = group_beams[g][y][1]
                y_sites = group_beams[g][y][2]

                merged_motifs = x_motifs + y_motifs

                merged_sites = []
                for s in range (len(x_sites)):
                    merged_sites.append(x_sites[s])
                merged_sites[a] = y_sites[a]
                merged_sites[b] = y_sites[b]

                score = get_score(build_ppm_from_brute_force(merged_motifs))
                joined.append([score, merged_motifs, merged_sites])

        beam = keep_top_n(joined, beam_width)

    # Sequential tail: extend every beam entry by every l-mer of each
    # remaining sequence, keep top beam_width (unchanged from topN).
    seeded_count = 2 * num_groups
    for tt in range(seeded_count, t):
        seq_index = order[tt]
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



