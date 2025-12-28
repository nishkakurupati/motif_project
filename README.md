# Motif Finding Project

## Overview

This project implements and evaluates a **motif-finding algorithm** using **synthetic benchmark data**.  
A *motif* is a short, recurring DNA pattern that serves as a biological signal for proteins interacting with DNA. Identifying motifs helps uncover regulatory regions and shared functional elements across sequences.

The goal of this project is to assess how accurately and efficiently a motif-finding algorithm can recover a **planted motif** under controlled experimental conditions.

The project is structured into **two levels**:
- **Level 1:** Perfect motifs (no noise)
- **Level 2:** Motifs with noise (1 character error)

This repository focuses on **Level 1**.

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

## Future Work

- Extend the project to **Level 2** with noisy motifs
- Explore optimized or probabilistic motif-finding algorithms
- Apply the method to real biological datasets

---
