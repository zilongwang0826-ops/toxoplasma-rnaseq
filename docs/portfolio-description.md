# Project description

Université Paris-Saclay, Master 2 workflow project | October 2026

This project reanalyses 12 public mouse-brain RNA-seq libraries and compares the results with Hu et al. (2020). A Snakemake workflow combines read-quality checks, Salmon quantification, tximport/DESeq2 analysis and gene-level GO enrichment. Comparisons use three infected and three control libraries at each of 11 and 33 days.

The analysis identifies 537 and 1,191 differential genes and recovers 74.36% and 78.04% of the paper's selected transcript lists. The main work is a documented analysis and comparison: gene and transcript results are kept separate, agreement is measured using transcript IDs, and a focused check examines the effect of two overlapping annotation records on GO counts. Source records, software versions, code and saved results are available in the repository.

These are whole-brain observations from an existing study. Pooling language leaves mouse-level independence unresolved, and selected-list recovery is not genome-wide accuracy. No new animal experiment or independent biological validation is claimed.
