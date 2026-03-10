# Motif Finding Project

## Overview

This project implements and evaluates a **motif-finding algorithm** using **synthetic benchmark data**.  
A *motif* is a short, recurring DNA pattern that serves as a biological signal for proteins interacting with DNA. Identifying motifs helps uncover regulatory regions and shared functional elements across sequences.

The goal of this project is to assess how accurately and efficiently a motif-finding algorithm can recover a **planted motif** under controlled experimental conditions.

The project is structured into **two levels**:
- **Level 1:** Perfect motifs (no noise)
- **Level 2:** Motifs with noise (some sequences missing the motif, others with a 1-character mutation)

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

