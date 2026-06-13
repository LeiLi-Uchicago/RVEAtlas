# RSVExplorer Methods and Information

## Note

While we included partial genomes to maximize G gene coverage, users should exercise caution when interpreting insertions or deletions (indels) in this region. Due to the hypervariable nature and repetitive sequences of the G gene, indels identified near the gene boundaries (N- and C-termini) may be alignment artifacts rather than true biological variants. We recommend manual verification of these gaps against the raw mapping data.

---

## Overview

RSVExplorer is an interactive Shiny application for exploring amino acid variation across Respiratory Syncytial Virus (RSV) genomes. It is designed to help users examine subtype-specific sequence diversity, compare clades, inspect variation at individual positions, and visualize gene-wide patterns such as conservation and fixed mutations.

The current version of the app is loaded with curated RSV-A and RSV-B amino acid usage data organized by gene, clade, year, and year-month. These data are derived from the raw count tables under `data/RSVA/output_tables` and `data/RSVB/output_tables`, with metadata loaded from the corresponding cleaned metadata tables for each subtype.

To improve performance on large usage tables, RSVExplorer uses a DuckDB-backed cache for interactive amino acid queries when DuckDB is available. The app still keeps the legacy RDS cache path as a fallback. Single Position Explorer queries are pushed down to DuckDB, including optional year-month filtering when viewing clade-grouped data.

The Genetic Clade tab provides a metadata-level view of whole-genome RSV clades. It uses a compact metadata summary cache to rank clade choices by sequence count, summarize active monthly periods and peak prevalence, plot monthly clade prevalence, and show top country, region, and host breakdowns without keeping additional full metadata copies in memory.

RSVExplorer is intended for researchers, bioinformaticians, genomic epidemiologists, and other users who want to investigate RSV evolutionary patterns, mutation dynamics, and lineage-specific amino acid changes through an accessible visual interface.

---

## Data Processing

We process nucleotide sequences by first cleaning and matching them to metadata, then splitting them into manageable chunks and running Nextclade3 on each chunk. From the Nextclade outputs, we parse per-protein amino acid states, codons, deletions, insertions, missing regions, and metadata groupings, then aggregate those parsed rows into protein-level usage tables by time and clade. Low-support insertions are tracked separately and excluded from the main usage tables unless they meet the configured support threshold.

Reference for RSV A: A/England/397/2017 (EPI_ISL_412866), GeneBank: PP109421
Reference for RSV B: B/Australia/VIC-RCH056/2019 (EPI_ISL_1653999), GeneBank: OP975389

---

## Update Log

### 2026-05-12

- Added Clade annotation and codon.

### 2026-05-05

- Added a top-level Genetic Clade tab for exploring whole-genome RSV clades. Added compact metadata-cache summaries for ranked clade choices, monthly prevalence, peak periods, and top metadata breakdowns.
- Updated usage-cache construction to read raw count tables from `RSVA/output_tables` and `RSVB/output_tables`.

### 2026-04-28

- Added DuckDB-backed usage query support to reduce memory pressure during interactive filtering, while preserving the RDS fallback path.
- Added optional year-month filtering for clade-grouped Single Position Explorer views.
- Updated Dataset Insights Custom Dataset Breakdown to require a specific subtype and use subtype-specific clade color palettes.

### 2026-04-10

- Updated data processing and UI to better present in-frame insertion/deletions.

### 2026-04-06

- Updated data using latest NextStrain data (Downloaded on 2026-04-06).

### 2026-03-06

- 1st version on-line. Using NextStrain data (Downloaded on 2025).
