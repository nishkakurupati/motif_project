# Create a list with all the different combinations
ml = [6, 7, 8]
sc = [5, 10, 20]
sl = [500, 1000, 2000]

# --- Gibbs + hard-EM engine knobs (tunable) ---
# Gibbs sweeps per restart (each sweep resamples every sequence's start once).
# Kept short: the sweeps explore, hard-EM does the polishing, so spending the
# budget on MORE restarts (more basins) beats longer individual runs.
gibbs_sweeps = 40
# Laplace pseudocount added to every base in every column when building a PPM.
pseudocount = 0.5
# Max hard-EM polish iterations after the Gibbs phase (stops early on convergence).
em_max_iters = 50


def restarts_for(sl):
    """More restarts for short sequences (cheap, and finding the true global
    optimum is the whole game); fewer for long sequences to bound runtime."""
    if sl <= 500:
        return 150
    if sl <= 1000:
        return 80
    return 40

params_list = []

params_list.append([6,10,500])
params_list.append([7,10,500])
params_list.append([8,10,500])
params_list.append([8,5,500])
params_list.append([8,20,500])
params_list.append([8,10,1000])
params_list.append([8,10,2000])
