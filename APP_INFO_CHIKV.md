# Amino Acid Explorer Methods and Information

---

## Overview

Amino Acid Explorer is an interactive Shiny application for exploring amino acid variation across viral genomes. It is designed to help users examine dataset-specific sequence diversity, compare clades or other groups, inspect variation at individual positions, and visualize gene-wide patterns such as conservation and fixed mutations.

The app reads a universal raw-data layout, which is the output of NextAA pipeline, under `data/raw`: metadata from `cleaned_metadata.tsv` and amino acid usage tables from `output_tables/usage_<gene>_by_<grouping>.tsv`. Genes and grouping factors are detected automatically from those filenames, so the same app can be rebuilt for a different virus species when the incoming files follow the same format.

To improve performance on large usage tables, the app builds both RDS caches and a DuckDB-backed cache for interactive amino acid queries when DuckDB is available. Single Position Explorer queries are pushed down to DuckDB, including optional year-month filtering when viewing clade-grouped data.

The Genetic Clade tab provides a metadata-level view of clades or compatible group annotations. It uses a compact metadata summary cache to rank clade choices by sequence count, summarize active monthly periods and peak prevalence, plot monthly clade prevalence, and show top country, region, and host breakdowns without keeping additional full metadata copies in memory.

Amino Acid Explorer is intended for researchers, bioinformaticians, genomic epidemiologists, and other users who want to investigate viral evolutionary patterns, mutation dynamics, and lineage-specific amino acid changes through an accessible visual interface.

---

## Data Processing

The upstream pipeline, NextAA, provides one cleaned metadata table and one or more count tables with at least `position`, `aa`, `codon`, and `count` columns. Group-specific files should also include the grouping column named in the filename, such as `clade` for `usage_C_by_clade.tsv`. Time-aware files should include `Year` and `Month`; the app derives `Year_Month` and `Year` views from those fields.

---

## Update Log

### 2026-05-06

- Converted the app to a universal raw-data loader. It discovers genes and grouping factors from `data/raw/output_tables`, builds RDS and DuckDB caches automatically.

