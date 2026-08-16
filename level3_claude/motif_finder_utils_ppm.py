import random
import math

# Motif finder for the OOPS model (One Occurrence Per Sequence): every sequence
# contains exactly one length-ml motif instance sampled from a shared PPM whose
# preferred base per column has probability motif_p (0.8 here), with random ACGT
# elsewhere. We recover per-sequence start positions in two stages:
#
#   SEARCH  -- Gibbs sampling + hard-EM on a log-odds objective, many restarts.
#              Log-odds is a strong search landscape (finds tight optima fast).
#   SELECT  -- among the diverse optima found across restarts, pick the one that
#              best fits a Binomial(n, motif_p) column profile. Raw log-odds
#              overfits to spurious perfectly-conserved columns (in 500bp x 10
#              seqs, chance near-matches can assemble a phantom motif tighter
#              than the real p=0.8 one); the Binomial profile penalizes that
#              over-conservation and prefers a genuine p=0.8-shaped alignment.

BASES = ['A', 'C', 'G', 'T']
BASE_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
LOG2_4 = 2.0  # log2(1 / 0.25) -- uniform background correction


def get_sequences_from_file(file_name):
    sequences = []

    filename = 'data_set/file_' + file_name + "sequences.fa"
    file_obj = open(filename, 'r')
    for line in file_obj:
        sequences.append(line.rstrip('\n'))
        sl = len(line.rstrip('\n'))
    file_obj.close()

    return sl, sequences


def build_ppm_from_brute_force(motifs):
    """Raw-frequency PPM from a list of equal-length motif strings.
    Kept for output/KL comparability with the other level-3 variants."""
    ppm = []
    for i in range(len(motifs[0])):
        dict_counter = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for motif in motifs:
            dict_counter[motif[i]] += 1
        for key in dict_counter:
            dict_counter[key] = dict_counter[key] / len(motifs)
        ppm.append(dict_counter)
    return ppm


def get_score(ppm):
    """Information content of a raw-frequency PPM (kept for compatibility)."""
    total_score = 0
    for column in ppm:
        for probability in column.values():
            if probability > 0:
                total_score += probability * math.log2(4 * probability)
    return total_score


# ----------------------------------------------------------------------------
# Engine internals (integer-encoded for speed)
# ----------------------------------------------------------------------------

def encode(sequence):
    """Map a DNA string to a list of base indices 0..3."""
    return [BASE_INDEX[c] for c in sequence]


def build_col_counts(starts, enc_seqs, ml):
    """Per-column base counts [ml][4] from current start positions."""
    counts = [[0, 0, 0, 0] for _ in range(ml)]
    for s in range(len(enc_seqs)):
        seq = enc_seqs[s]
        p = starts[s]
        for i in range(ml):
            counts[i][seq[p + i]] += 1
    return counts


def _leave_one_out(counts, seq, p_old, ml):
    for i in range(ml):
        counts[i][seq[p_old + i]] -= 1


def _add_window(counts, seq, p_new, ml):
    for i in range(ml):
        counts[i][seq[p_new + i]] += 1


# ---- log-odds scoring (SEARCH objective) -----------------------------------

def counts_to_logodds(counts, n, ml, alpha):
    """Log-odds weights [ml][4]: log2((count+alpha)/(n+4alpha) / 0.25)."""
    denom = n + 4 * alpha
    lw = [[0.0, 0.0, 0.0, 0.0] for _ in range(ml)]
    for i in range(ml):
        for b in range(4):
            lw[i][b] = math.log2((counts[i][b] + alpha) / denom) + LOG2_4
    return lw


def _window_logodds(seq, lw, ml, num_positions):
    scores = []
    for p in range(num_positions):
        total = 0.0
        for i in range(ml):
            total += lw[i][seq[p + i]]
        scores.append(total)
    return scores


def logodds_alignment_score(starts, enc_seqs, ml, alpha):
    counts = build_col_counts(starts, enc_seqs, ml)
    lw = counts_to_logodds(counts, len(enc_seqs), ml, alpha)
    total = 0.0
    for s in range(len(enc_seqs)):
        seq = enc_seqs[s]
        p = starts[s]
        for i in range(ml):
            total += lw[i][seq[p + i]]
    return total


# ---- Binomial-profile scoring (SELECTION discriminator) --------------------

def binomial_logpmf_table(n, p):
    """logpmf[k] = log Binomial(k; n, p) for k = 0..n. Peaks near k = n*p and
    penalizes over-conservation, so it favours genuine p-shaped columns."""
    logp = math.log(p)
    logq = math.log(1.0 - p)
    table = []
    for k in range(n + 1):
        logc = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
        table.append(logc + k * logp + (n - k) * logq)
    return table


def profile_score(starts, enc_seqs, ml, logpmf):
    counts = build_col_counts(starts, enc_seqs, ml)
    total = 0.0
    for i in range(ml):
        total += logpmf[max(counts[i])]
    return total


# ---- search moves ----------------------------------------------------------

def sample_from_scores(scores, beta):
    """Sample an index proportional to 2^(beta*(score - max))."""
    max_score = max(scores)
    weights = []
    total = 0.0
    for s in scores:
        w = 2.0 ** (beta * (s - max_score))
        weights.append(w)
        total += w
    r = random.random() * total
    cum = 0.0
    for i in range(len(weights)):
        cum += weights[i]
        if r <= cum:
            return i
    return len(weights) - 1


def hard_em(starts, enc_seqs, ml, alpha, num_positions, max_iters):
    """Iterate {rebuild log-odds from other seqs -> set each start = argmax
    window} with leave-one-out until starts stop changing."""
    sc = len(enc_seqs)
    starts = list(starts)
    counts = build_col_counts(starts, enc_seqs, ml)
    for _ in range(max_iters):
        changed = False
        for s in range(sc):
            seq = enc_seqs[s]
            p_old = starts[s]
            _leave_one_out(counts, seq, p_old, ml)
            lw = counts_to_logodds(counts, sc - 1, ml, alpha)
            scores = _window_logodds(seq, lw, ml, num_positions)
            best_p = 0
            best_v = scores[0]
            for p in range(1, num_positions):
                if scores[p] > best_v:
                    best_v = scores[p]
                    best_p = p
            starts[s] = best_p
            _add_window(counts, seq, best_p, ml)
            if best_p != p_old:
                changed = True
        if not changed:
            break
    return starts


def gibbs_restart(enc_seqs, ml, alpha, num_positions, sweeps):
    """One Gibbs run: random init, resample each start from the leave-one-out
    log-odds distribution, temperature cooled (beta rising) across sweeps."""
    sc = len(enc_seqs)
    starts = [random.randint(0, num_positions - 1) for _ in range(sc)]
    counts = build_col_counts(starts, enc_seqs, ml)
    for sweep in range(sweeps):
        beta = 0.5 + 1.5 * (sweep / (sweeps - 1)) if sweeps > 1 else 1.0
        for s in range(sc):
            seq = enc_seqs[s]
            p_old = starts[s]
            _leave_one_out(counts, seq, p_old, ml)
            lw = counts_to_logodds(counts, sc - 1, ml, alpha)
            scores = _window_logodds(seq, lw, ml, num_positions)
            p_new = sample_from_scores(scores, beta)
            starts[s] = p_new
            _add_window(counts, seq, p_new, ml)
    return starts


def find_motif_sites(sequences, sl, ml, num_restarts, gibbs_sweeps,
                     pseudocount, em_max_iters, motif_p=0.8, select='profile'):
    """Recover per-sequence motif start positions. Returns (best_ppm, sites)
    with sites in ORIGINAL sequence order.

    Search with log-odds (many Gibbs+EM restarts), collect the distinct optima,
    then choose the final alignment with `select`:
      'profile'  -> Binomial(n, motif_p) column profile (default; anti-overfit)
      'logodds'  -> raw log-odds (original maximum-likelihood pick)
    """
    enc_seqs = [encode(seq) for seq in sequences]
    num_positions = sl - ml + 1
    n = len(enc_seqs)
    logpmf = binomial_logpmf_table(n, motif_p)

    seen = set()
    candidates = []
    for _ in range(num_restarts):
        starts = gibbs_restart(enc_seqs, ml, pseudocount, num_positions, gibbs_sweeps)
        starts = hard_em(starts, enc_seqs, ml, pseudocount, num_positions, em_max_iters)
        key = tuple(starts)
        if key not in seen:
            seen.add(key)
            candidates.append(starts)

    if select == 'logodds':
        scorer = lambda st: logodds_alignment_score(st, enc_seqs, ml, pseudocount)
    else:
        scorer = lambda st: profile_score(st, enc_seqs, ml, logpmf)

    best_starts = candidates[0]
    best_val = scorer(best_starts)
    for st in candidates[1:]:
        v = scorer(st)
        if v > best_val:
            best_val = v
            best_starts = st

    motifs = [sequences[s][best_starts[s]:best_starts[s] + ml] for s in range(len(sequences))]
    best_ppm = build_ppm_from_brute_force(motifs)
    return best_ppm, best_starts


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
    for i in range(sc):
        file_obj.write(str(sites[i]) + "\n")
    file_obj.close()
