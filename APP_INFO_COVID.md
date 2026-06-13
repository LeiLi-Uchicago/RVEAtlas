# COVIDExplorer Methods and Information

---

## Overview

COVIDExplorer is an interactive Shiny application for exploring amino acid variation across SARS-COVID-2 genomes. It is designed to help users examine subtype-specific sequence diversity, compare clades, inspect variation at individual positions, and visualize gene-wide patterns such as conservation and fixed mutations.

The current version of the app is loaded with curated COVID amino acid usage data organized by gene, clade, year, and year-month. 

COVIDExplorer is intended for researchers, bioinformaticians, genomic epidemiologists, and other users who want to investigate RSV evolutionary patterns, mutation dynamics, and lineage-specific amino acid changes through an accessible visual interface.

---

## Data Processing

We process nucleotide sequences by first cleaning and matching them to metadata, then splitting them into manageable chunks and running Nextclade3 on each chunk. From the Nextclade outputs, we parse per-protein amino acid states, codons, deletions, insertions, missing regions, and metadata groupings, then aggregate those parsed rows into protein-level usage tables by time and clade. Low-support insertions are tracked separately and excluded from the main usage tables unless they meet the configured support threshold.

Reference for COVID: Wuhan-Hu-1/2019, GeneBank: MN908947

---

## Update Log

### 2026-05-22

- Update clade and codon.

### 2026-05-05

- Added the Genetic Clade tab. Users can select a clade annotation, search or choose a ranked clade, and review clade-level summary statistics.

### 2026-04-28

- 1st version on-line. Using latest NextStrain data (Downloaded on 2026-04-06).
