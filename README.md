# Motif Finding Project

## Overview

This project implements and evaluates a **motif-finding algorithm** using **synthetic benchmark data**.  
A *motif* is a short, recurring DNA pattern that serves as a biological signal for proteins interacting with DNA. Identifying motifs helps uncover regulatory regions and shared functional elements across sequences.

The goal of this project is to assess how accurately and efficiently a motif-finding algorithm can recover a **planted motif** under controlled experimental conditions.

The project is structured into **three levels** of increasing difficulty:
- **Level 1:** Perfect motifs (no noise) — the exact motif appears in every sequence
- **Level 2:** Motifs with noise (some sequences missing the motif, others with a 1-character mutation)
- **Level 3:** Probabilistic motifs (PPM) — every sequence receives an independently sampled instance, so no exact motif string exists

Level 3 is explored through **five algorithm variants**, each in its own directory:

| Directory | Approach |
|-----------|----------|
| `level3_simple_1pass` | Deterministic greedy, single pass |
| `level3_randomized` | Greedy with 50 shuffled-order restarts |
| `level3_beam_search_2way` | Beam search (width 100), 2-sequence seed |
| `level3_beam_search_6way` | Beam search with 3-pair seeding and cross-product join |
| `level3_claude` | Gibbs sampling + hard-EM, with Binomial-profile selection |

---

## What Is a Motif?

A motif is a short nucleotide pattern (typically a few base pairs long) that appears repeatedly within or across DNA sequences.  
Despite their small size, motifs play a crucial role in gene regulation by providing recognizable binding sites for proteins.

Motifs also act as **signatures** in DNA sequences, making them useful for:
- Identifying shared functional regions
- Comparing multiple sequences
- Discovering regulatory elements

---

## Project Description

The project consists of three main components:

1. **Benchmark Construction**
2. **Motif Finder Implementation**
3. **Evaluation and Performance Analysis**

---

## 1. Benchmark Construction

Synthetic benchmark datasets are generated to enable controlled and reproducible evaluation.

Each dataset is defined by:
- **ML** – Motif Length  
- **SL** – Sequence Length  
- **SC** – Number of Sequences  

### Dataset Generation
For each dataset:
- `SC` random DNA sequences of length `SL` are generated with uniform nucleotide frequencies
- A random motif of length `ML` is created
- The motif is planted at a random position in each sequence

Each dataset directory contains:
- `sequences.fa` – DNA sequences (FASTA format)
- `sites.txt` – True motif locations
- `motif.txt` – Planted motif
- `motiflength.txt` – Motif length

---

## Experimental Design

### Default Parameters
- ML = 8  
- SL = 500  
- SC = 10  

### Parameter Variations
To study the effect of individual parameters, one parameter is varied at a time:

| Parameter | Values |
|---------|--------|
| Motif Length (ML) | 6, 7, 8 |
| Sequence Count (SC) | 5, 10, 20 |
| Sequence Length (SL) | 500, 1000, 2000 |

- 10 independent datasets are generated for each configuration  
- Total number of datasets: **70**

---

## 2. Motif Finder

The motif-finding program uses only:
- `sequences.fa`
- `motiflength.txt`

### Output
For each dataset, the program produces:
- `predictedmotif.txt` – Predicted motif sequence
- `predictedsites.txt` – Predicted motif locations

Both files follow the same format as the corresponding ground-truth files.

The algorithm uses a **brute-force approach**, iteratively eliminating candidate motifs that do not appear in all sequences. Because the benchmark contains **perfect motifs**, an exact match is guaranteed.

---

## 3. Evaluation and Analysis

The evaluation script runs the motif finder on all benchmark datasets and measures:

- **Motif Accuracy** – Exact match between true and predicted motif
- **Site Accuracy** – Overlap between predicted and true motif locations
- **Runtime** – Execution time in milliseconds

Results are averaged across the 10 datasets for each parameter configuration and visualized using plots.

---

## Results

### Varying Motif Length (ML)
*(SL = 500, SC = 10)*

| ML | Runtime (ms) | Motif Accuracy (%) | Site Accuracy (%) |
|----|-------------|-------------------|------------------|
| 6  | 47.83 | 100 | 96 |
| 7  | 46.32 | 100 | 99 |
| 8  | 42.72 | 100 | 100 |

---

### Varying Sequence Length (SL)
*(ML = 8, SC = 10)*

| SL | Runtime (ms) | Motif Accuracy (%) | Site Accuracy (%) |
|----|-------------|-------------------|------------------|
| 500  | 43.31 | 100 | 100 |
| 1000 | 183.38 | 100 | 99 |
| 2000 | 758.73 | 100 | 99 |

---

### Varying Sequence Count (SC)
*(SL = 500, ML = 8)*

| SC | Runtime (ms) | Motif Accuracy (%) | Site Accuracy (%) |
|----|-------------|-------------------|------------------|
| 5  | 42.62 | 100 | 100 |
| 10 | 43.31 | 100 | 100 |
| 20 | 44.27 | 100 | 99.5 |

---

## Discussion

### Motif Matching Performance
Motif matching accuracy was **100%** across all parameter combinations. This is expected because the motifs were analyzed without noise, ensuring exact matches.

### Site Matching Performance
- **Sequence Length:** Longer sequences increase the likelihood of duplicate motif occurrences, reducing measured site accuracy.
- **Motif Length:** Shorter motifs are more ambiguous and appear more frequently by chance, lowering site accuracy.
- **Sequence Count:** Increasing the number of sequences had minimal impact on site accuracy.

### Runtime Performance
- **Sequence Length:** Had the largest impact on runtime due to the increased number of possible motif candidates.
- **Sequence Count:** Runtime increased only slightly, as candidate motifs converge quickly.
- **Motif Length:** Small variations had limited impact on runtime.

---

## Conclusion

This project demonstrates that a simple brute-force motif-finding algorithm can achieve **perfect motif recovery** under ideal, noise-free conditions. However, runtime performance degrades with longer sequences, highlighting the need for more scalable approaches when handling large or noisy datasets.

---

---

# Level 2: Motif Finder with Noise

## Overview

Level 2 extends the benchmark to a more realistic, noisy setting. The planted motif is no longer guaranteed to appear perfectly in every sequence — some sequences may be missing it entirely, and others may contain a mutated copy with one character changed.

---

## 1. Benchmark Construction

Dataset generation follows the same structure as Level 1 (same parameter combinations, 10 repetitions each, 70 datasets total), with one key difference in how the motif is planted:

### Noise Model

For each sequence, the motif is inserted as follows:

- **~20% chance:** The motif is **not inserted** at all (position recorded as `-1`)
- **~80% chance:** The motif is inserted, but with a further split:
  - **~20% of insertions:** A **single random character** in the motif is replaced with a different nucleotide (1-character mutation)
  - **~80% of insertions:** The **exact motif** is inserted

This means across a dataset, the planted motif appears exactly in the majority of sequences, with some sequences holding a 1-character variant and some holding no motif at all.

Each dataset directory contains the same four files as Level 1:
- `sequences.fa` – DNA sequences (FASTA format)
- `sites.txt` – True motif locations (`-1` where no motif was planted)
- `motif.txt` – Planted motif
- `motiflength.txt` – Motif length

---

## 2. Motif Finder

Because the exact motif is no longer guaranteed to appear in every sequence, the Level 1 brute-force elimination approach breaks down. Level 2 uses a **frequency-based (majority vote) algorithm** instead.

### Algorithm

1. For each sequence, extract all substrings of length `ML` and store them in a per-sequence dictionary (only the first occurrence of each substring is recorded, along with its position)
2. Count how many sequences each unique substring appears in
3. Select the substring with the highest count as the predicted motif
4. Look up the predicted motif's position in each sequence's dictionary; report `-1` if it does not appear

This approach works because the exact motif appears in the majority of sequences, making it the most frequent length-`ML` substring across the dataset.

### Output

For each dataset, the program produces:
- `predictedmotif.txt` – Predicted motif sequence
- `predictedsites.txt` – Predicted motif locations (`-1` where the motif was not found)

---

## 3. Evaluation and Analysis

The evaluation script (`evaluate_with_noise.py`) runs the motif finder on all 70 benchmark datasets and measures the same three metrics as Level 1:

- **Motif Accuracy** – Exact match between true and predicted motif
- **Site Accuracy** – Overlap between predicted and true motif locations
- **Runtime** – Execution time in milliseconds

Results are averaged across the 10 datasets for each parameter configuration.

---

## Results

### Varying Motif Length (ML)
*(SL = 500, SC = 10)*

| ML | Runtime (ms) | Motif Accuracy (%) | Site Accuracy (%) |
|----|-------------|-------------------|------------------|
| 6  | 0.00171387 | 100| 90 |
| 7  | 0.00168307 | 90 | 85 |
| 8  | 0.00166142 | 100 | 95 |

---

### Varying Sequence Length (SL)
*(ML = 8, SC = 10)*

| SL | Runtime (ms) | Motif Accuracy (%) | Site Accuracy (%) |
|----|-------------|-------------------|------------------|
| 500  | 0.00166142 |100 | 95 |
| 1000 | 0.00324888 | 100 | 92 |
| 2000 | 0.00637912 | 100| 93 |

---

### Varying Sequence Count (SC)
*(SL = 500, ML = 8)*

| SC | Runtime (ms) | Motif Accuracy (%) | Site Accuracy (%) |
|----|-------------|-------------------|------------------|
| 5  | 0.00093607 | 90 | 78|
| 10 | 0.00166142 | 100 | 95 |
| 20 | 0.00316435| 100 | 90.5 |


---

## Discussion

### Motif Matching Performance
Unlike Level 1, motif accuracy is no longer guaranteed to be 100%. The majority-vote algorithm recovers the correct motif only when the planted motif is the most frequent length-`ML` substring across the dataset. With shorter motifs (smaller `ML`), random substrings are more likely to appear frequently by chance, increasing the risk of a false winner. With fewer sequences (`SC`), the signal from the planted motif is weaker relative to noise.

### Site Matching Performance
Site accuracy is expected to be lower than Level 1 for two reasons:
1. Some sequences had no motif planted (ground truth = `-1`), so those sites cannot be matched
2. Some sequences received a mutated copy, so the exact predicted motif will not be found there, resulting in a predicted position of `-1`

### Runtime Performance
The frequency-based algorithm scans each sequence once to build its substring dictionary, making it more efficient than the Level 1 elimination approach for longer sequences. Runtime is expected to scale more gradually with increasing `SL`.

---

## Conclusion

Level 2 demonstrates that motif finding under realistic noise conditions is significantly harder. The majority-vote approach provides a practical heuristic that leverages the signal from the majority of sequences, but it is sensitive to the balance between motif length, sequence count, and noise level.

---

---

# Level 3: Probabilistic Motifs (PPM)

## Overview

Levels 1 and 2 both plant a single motif **string**, so there is always one exact answer to recover. Level 3 replaces that with a **Position Probability Matrix (PPM)**: the motif is a probability distribution over `A/C/G/T` for each of its `ML` columns, and every sequence receives an **independently sampled instance** drawn from that matrix.

This is the realistic formulation of the problem, and it is fundamentally harder:

- No two planted instances are necessarily identical, so there is no exact string to search for
- The exact-match and majority-vote algorithms from Levels 1 and 2 do not apply at all
- The planted alignment is no longer guaranteed to be the *highest scoring* alignment — random windows can assemble into a tighter-looking motif than the real one

Because of this, Level 3 is explored through **five algorithm variants**, each in its own self-contained directory. They share identical benchmark, scoring, and evaluation code; only the finder changes.

---

## 1. Benchmark Construction

The parameter grid is unchanged from Levels 1 and 2 (same 7 configurations, 10 repetitions each, 70 datasets total). What changes is how the motif is defined and planted.

### PPM Generation

For each column of the motif, one nucleotide is chosen at random as the **preferred base** and given probability `p`; the remaining three each get `(1 - p) / 3`:

| Base | Probability (p = 0.8) |
|------|----------------------|
| Preferred | 0.8 |
| Each other | 0.0667 |

All results below use **p = 0.8**.

### Planting

For each sequence, a fresh motif instance is sampled column by column from the PPM and inserted at a random position. With `ML = 8` and `p = 0.8`, only about `0.8^8 ≈ 17%` of instances are the full consensus string — the rest each differ from it in one or more positions.

Unlike Level 2, **every sequence receives an instance** (the One Occurrence Per Sequence, or OOPS, model). There are no `-1` sites; all the difficulty comes from the variability of the instances, not from missing ones.

### Files

Same four files per dataset, but `motif.txt` now holds the PPM as a 4-row matrix rather than a string:

```
A 0.8 0.0667 0.0667 ...
T 0.0667 0.0667 0.8 ...
C 0.0667 0.8 0.0667 ...
G 0.0667 0.0667 0.0667 ...
```

The finder writes `predicted_motif.txt` (same matrix format) and `predicted_sites.txt`.

> Note the underscore: Level 3 outputs are `predicted_motif.txt` / `predicted_sites.txt`, whereas Levels 1 and 2 use `predictedmotif.txt` / `predictedsites.txt`.

---

## 2. Evaluation Metrics

Because the answer is now a matrix, exact motif matching is meaningless and is replaced by a divergence measure.

### Motif Accuracy → KL Divergence

The predicted PPM is compared to the true PPM using **Kullback-Leibler divergence**, summed over every column and base. **Lower is better**; 0 would mean a perfect match.

The predicted PPM is built as a raw frequency count over the `SC` recovered instances, so with `SC = 10` its entries are multiples of 0.1 and many are exactly 0, while the true PPM has entries of 0.8 and 0.0667. This representational mismatch means even a *perfectly aligned* prediction carries a substantial KL penalty. KL values are therefore useful for **comparing variants against each other**, not as an absolute measure of correctness.

### Site Accuracy

Unchanged in spirit from Levels 1 and 2: the fraction of sequences where the predicted start position exactly equals the true one. This is the **primary metric** for Level 3 — it is the only one that directly answers "did we find the real motif?"

### Runtime

Seconds per dataset. Note the unit change: Level 3 runtimes are **seconds**, not milliseconds.

---

## 3. Motif Finder Variants

All five variants score a candidate alignment with the same objective — the **information content** of the resulting PPM (`get_score`) — except where noted. They differ in how they search the space of alignments.

### 3.1 `level3_simple_1pass`

A deterministic greedy search, and the baseline everything else is measured against.

1. Exhaustively try **all pairs** of l-mers from sequences 0 and 1, keep the best-scoring pair
2. For each remaining sequence in order, add the single l-mer that most improves the score
3. Never revisit an earlier choice

One pass, no randomness. Fast, but every later decision is locked to the initial pair, so a bad seed cannot be recovered from.

### 3.2 `level3_randomized`

The same greedy algorithm, wrapped in **random restarts**.

- The sequence order is **shuffled** before each trial, so a different pair seeds the search each time
- `num_trials = 50` trials per dataset; the highest-scoring trial wins
- Sites are tracked by **original sequence index** throughout, so shuffling does not corrupt the reported positions

This addresses the "bad seed" failure of the single pass by sampling 50 different starting points.

### 3.3 `level3_beam_search_2way`

Replaces greedy's "keep the single best" with **beam search**: keep the top `N` partial alignments at every step.

- `beam_width = 100` — at each extension, every surviving partial alignment is extended by every candidate l-mer, and the top 100 results survive
- `num_trials = 5`, reduced from 50 to offset the much higher cost per trial
- Seeded from all pairs of the first **2** sequences in the shuffled order (hence *2way*)

Greedy is the special case `beam_width = 1`. Widening the beam lets a partial alignment that looks mediocre early survive long enough to prove itself.

### 3.4 `level3_beam_search_6way`

Beam search with **multi-group seeding**, to reduce dependence on which two sequences happen to seed the search.

- The shuffled order is split into up to **3 disjoint pairs** (up to 6 sequences, hence *6way*)
- Each pair is seeded exhaustively and keeps its **own** top-100 beam
- The three beams are then **cross-product joined** and pruned back to 100, so the seed reflects evidence from 6 sequences rather than 2
- Remaining sequences are extended sequentially exactly as in the 2way variant
- Falls back to single-sequence seeding when there are too few sequences to form a pair

### 3.5 `level3_claude`

A rewrite on a different algorithmic foundation: **Gibbs sampling with hard-EM refinement**, the approach used by real motif finders such as MEME and the Gibbs Motif Sampler.

**Search.** Each restart begins from random start positions and then:
- **Gibbs phase** — repeatedly resample each sequence's start position from the leave-one-out **log-odds** distribution over windows, with the temperature cooled across sweeps (`gibbs_sweeps = 40`)
- **Hard-EM phase** — deterministically snap each start to its argmax window until nothing moves (`em_max_iters = 50`)
- Restarts are adaptive by sequence length: **150** for `SL ≤ 500`, **80** for `SL ≤ 1000`, **40** beyond

**Selection.** This is the variant's key idea. Searching on log-odds finds tight alignments, but the tightest alignment is often *not* the real one — in 500 bp × 10 sequences, chance near-matches can assemble a phantom motif more conserved than a genuine p = 0.8 one. So the final answer is chosen by a **different** objective from the one used to search: among the distinct optima collected across all restarts, pick the one whose column conservation best fits a **Binomial(SC, 0.8)** profile. That distribution peaks near `SC × 0.8` and penalizes *over*-conservation, which is exactly the failure mode raw log-odds exhibits.

**Implementation.** Sequences are integer-encoded, column counts are updated incrementally via leave-one-out add/remove rather than rebuilt, and PPMs use a Laplace `pseudocount = 0.5`. These make the many-restart budget affordable.

The directory also includes `run_subset.py`, a fast-iteration harness that runs a small `k`-range across all 7 configurations and reports site accuracy against a hardcoded `INCUMBENT` table of the best results from the earlier variants.

---

## 4. Results

All numbers are averaged over the 10 datasets per configuration, at **p = 0.8**, and are taken from each variant's `results_08.txt`. Because all five directories generate their benchmark with the same code and the same seed, they are evaluated on **identical datasets**.

The **default** column is the shared baseline configuration `ML = 8, SC = 10, SL = 500`.

### Site Accuracy (fraction of sites exactly matched — higher is better)

| Variant | ML=6 | ML=7 | default | SL=1000 | SL=2000 | SC=5 | SC=20 |
|---------|-----:|-----:|--------:|--------:|--------:|-----:|-----:|
| `simple_1pass`     | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| `randomized`       | 0.00 | 0.00 | 0.08 | **0.06** | 0.00 | 0.00 | 0.21 |
| `beam_search_2way` | 0.05 | **0.19** | 0.19 | **0.06** | 0.00 | 0.08 | 0.24 |
| `beam_search_6way` | 0.05 | 0.18 | 0.18 | 0.00 | 0.00 | 0.08 | 0.31 |
| `claude`           | 0.05 | 0.14 | **0.20** | 0.00 | 0.00 | 0.08 | **0.49** |

### Mean KL Divergence (lower is better)

| Variant | ML=6 | ML=7 | default | SL=1000 | SL=2000 | SC=5 | SC=20 |
|---------|-----:|-----:|--------:|--------:|--------:|-----:|-----:|
| `simple_1pass`     | 15.90 | 17.97 | 18.46 | 18.99 | 19.27 | 21.65 | 17.63 |
| `randomized`       | 17.27 | 17.39 | 17.28 | **15.26** | 20.68 | 20.94 | 10.77 |
| `beam_search_2way` | 15.69 | **10.28** | 14.52 | 18.39 | 19.63 | 20.21 | 10.72 |
| `beam_search_6way` | 15.70 | 12.13 | 15.11 | 19.83 | **18.31** | **18.62** | 7.64 |
| `claude`           | **15.30** | 11.83 | **14.26** | 19.68 | 20.73 | 18.86 | **4.15** |

### Runtime (seconds per dataset — lower is better)

| Variant | ML=6 | ML=7 | default | SL=1000 | SL=2000 | SC=5 | SC=20 |
|---------|-----:|-----:|--------:|--------:|--------:|-----:|-----:|
| `simple_1pass`     | **0.88** | **1.01** | **1.14** | **4.57** | **18.48** | **1.13** | **1.18** |
| `randomized`       | 48.60 | 55.19 | 62.17 | 257.25 | 1015.14 | 60.83 | 66.30 |
| `beam_search_2way` | 17.59 | 20.46 | 22.60 | 63.07 | 188.77 | 12.68 | 49.97 |
| `beam_search_6way` | 27.35 | 29.21 | 32.60 | 114.44 | 425.71 | 17.60 | 60.82 |
| `claude`           | 8.81 | 9.20 | 10.51 | 11.46 | 11.79 | 5.49 | 20.44 |

---

## Discussion

### Site Accuracy Is Low Everywhere

The headline result is that **no variant exceeds 0.49, and most sit below 0.20**. This is the honest outcome of the PPM formulation: with instances sampled at p = 0.8, the planted alignment frequently is not the best-scoring alignment available, so even a search that finds the global optimum of the objective will report the wrong sites. The ceiling here is set by the **objective function**, not by search quality.

This is the single biggest difference from Levels 1 and 2, where the planted answer was also the optimal one and the only question was whether the algorithm could find it.

### Search Quality Still Matters

The progression from `simple_1pass` to the later variants is nonetheless real and monotone at the default configuration:

| Variant | Site accuracy (default) |
|---------|------------------------|
| `simple_1pass` | 0.00 |
| `randomized` | 0.08 |
| `beam_search_2way` | 0.19 |
| `beam_search_6way` | 0.18 |
| `claude` | 0.20 |

`simple_1pass` scores **exactly zero everywhere** — a single deterministic greedy pass never recovers a correct site in any of the 70 datasets. Random restarts alone lift this to 0.08, and widening the beam roughly doubles it again. Beyond that the returns flatten: `6way` seeding does not beat `2way` at the default, and `claude` improves on it only slightly.

### Where the Variants Genuinely Separate: SC = 20

The clearest separation is at **SC = 20**, where site accuracy climbs 0.21 → 0.24 → 0.31 → **0.49** and KL falls 10.77 → 10.72 → 7.64 → **4.15**. More sequences mean more evidence per column, which both sharpens the true signal and makes it harder for a phantom alignment to look conserved by chance. This is precisely the regime `claude`'s Binomial-profile selection is designed for — the Binomial(SC, p) distribution becomes much more discriminating as SC grows — and it is where that idea pays off, nearly doubling the next-best variant.

Conversely, at **SC = 5** every variant collapses to 0.08 or below. Five samples are simply not enough evidence to distinguish a real p = 0.8 motif from noise, regardless of algorithm.

### Long Sequences Defeat Every Variant

At **SL = 2000, every variant scores 0.00**, and at SL = 1000 only `randomized` (0.06) and `beam_search_2way` (0.06) register anything at all. Each additional position is another chance for a spurious window to beat the planted one, and the number of candidate alignments grows as `SL^SC`. Interestingly the *more sophisticated* variants do worse here: a better search finds a better-scoring alignment, and at SL = 2000 the better-scoring alignment is reliably the wrong one. Better optimization of the wrong objective produces worse answers.

### Shorter Motifs Are Harder

`ML = 6` sits at 0.05 or below for every variant, against 0.14–0.20 at `ML = 8`. A 6-column motif carries less information and is more easily imitated by chance — the same effect observed in Levels 1 and 2, but far more pronounced once instances are probabilistic.

### Runtime

The exhaustive all-pairs seed used by the greedy and beam variants is **O(SL²)**, which dominates everything else. This is why `randomized` — 50 trials, each paying that cost — reaches **1015 seconds per dataset** at SL = 2000. `beam_search_2way` is *faster* than `randomized` despite far more work per trial, purely because it runs 5 trials instead of 50. `6way` costs roughly double `2way`, since it pays the all-pairs seed three times.

`claude` breaks this pattern. Gibbs plus hard-EM is **O(restarts × sweeps × SC × SL)** — linear in `SL`, with no all-pairs seed anywhere — so its runtime is nearly **flat** across sequence length (10.51 → 11.46 → 11.79 s as SL goes 500 → 1000 → 2000), helped further by the adaptive restart schedule. At SL = 2000 it is roughly **86× faster** than `randomized` and **36× faster** than `6way`, while matching or beating both on accuracy. Only `simple_1pass` is cheaper, and it never finds anything.

---

## Conclusion

Level 3 shows that the probabilistic (PPM) formulation is a categorically harder problem than exact or lightly-mutated motif planting. Progressively stronger search — random restarts, then beam search, then Gibbs sampling with hard-EM — produces steady gains, and the split-objective idea in `level3_claude` (search on log-odds, select on a Binomial column profile) gives the strongest results at high sequence counts while being dramatically cheaper on long sequences.

But the dominant limitation is not search. At SL ≥ 1000 and at SC = 5, every variant fails regardless of how thoroughly it explores, because the planted alignment is not the one the scoring function prefers. Meaningful further progress would require a better **objective** — a background model, a stronger prior on column composition, or explicit modelling of the p = 0.8 generative process — rather than a better optimizer.

---

## Running the Code

Each level directory is self-contained and uses **relative paths**, so scripts must be run from inside their own directory:

```bash
cd level3_claude
python3 create_dataset_ppm.py   # generate the 70 benchmark datasets into data_set/
python3 evaluate_ppm.py         # run the finder on all datasets and print all three metrics
```

For Levels 1 and 2 the equivalent entry points are `create_dataset.py` / `evaluate.py` and `create_dataset_noise.py` / `evaluate_with_noise.py`.

Dataset generation is seeded (`random.seed(0)`), and `level3_claude` additionally seeds its randomized search (`random.seed(1)`), so both are reproducible run to run. The other Level 3 variants do not seed their search, so their results will vary slightly between runs.

For fast iteration on `level3_claude` without running the full 70-dataset sweep:

```bash
cd level3_claude
python3 run_subset.py 3          # all 7 configurations, k = 0..2
python3 run_subset.py 5 logodds  # k = 0..4, using log-odds selection instead of the Binomial profile
```

