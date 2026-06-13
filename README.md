# FLU Amino Acid Divergence Explorer

![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)
![Shiny](https://img.shields.io/badge/Built_with-R_Shiny-success)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-purple)
![License](https://img.shields.io/badge/License-MIT-green)

**FLU Amino Acid Divergence Explorer** is an interactive Shiny application for exploring influenza amino acid variation across subtypes, clades, countries/regions, and collection dates. It supports dataset-level summaries, clade prevalence tracking, single-position residue distributions, pairwise clade comparisons, and gene-wide conservation/mutation views.

## Key Features

### Dataset Insights

- Summarizes total sequences, represented countries, and collection time span.
- Visualizes sequencing volume over time, regional composition, and subtype-specific metadata breakdowns.

### Genetic Clade

- Lets users choose a clade annotation and then search/select ranked clades.
- Shows total sequence count, rank/share, active period, peak prevalence month, monthly prevalence curves, and metadata breakdowns.

### Single Position Explorer

- Shows amino acid distributions at a selected gene position.
- Supports grouping by year, year-month, clade, and available metadata annotations.
- Can display percentages or raw counts and export plots/tables.

### Pairwise Comparison

- Compares two selected groups/clades across genes.
- Identifies fixed or near-fixed amino acid differences using a user-defined dominant-frequency threshold.

### Gene-Wide Landscapes

- **Conservation (Entropy):** Maps positional Shannon entropy across a gene.
- **Mutation Tracker (Lollipop):** Visualizes fixed amino acid mutations in a target group compared with a reference group.

## Local Setup

### 1. Install R

Install R 4.0 or newer from [CRAN](https://cran.r-project.org/). RStudio is optional but recommended for interactive development.

### 2. Clone the Repository

```bash
git clone https://github.com/LeiLi-Uchicago/FLUExplorer.git
cd FLUExplorer
```

### 3. Install R Packages

Open R from the project folder and install the required packages:

```r
install.packages(c(
  "shiny",
  "dplyr",
  "ggplot2",
  "DT",
  "readr",
  "tidyr",
  "openxlsx",
  "plotly",
  "waiter",
  "lubridate",
  "tidyverse",
  "shinyWidgets",
  "shinyjs",
  "viridis",
  "scales",
  "ggtext",
  "DBI",
  "duckdb"
))
```

`duckdb` and `DBI` are strongly recommended because the app uses a compact DuckDB usage cache to avoid loading large count tables into memory. If `duckdb` is not installed, the app falls back to legacy RDS lazy loading when those cache files are available.

### 4. Organize Raw Data

Place raw FLU metadata and count tables under `data/raw/FLU`. Each subtype should have its own folder:

```text
data/
└── raw/
    └── FLU/
        ├── H1N1/
        │   ├── metadata_merged_annotated.csv
        │   └── count/
        │       ├── HA/
        │       │   ├── aa_usage_by_HA_clade.csv
        │       │   ├── aa_usage_by_NA_clade.csv
        │       │   └── aa_usage_by_Year_Month.csv
        │       └── ...
        ├── H3N2/
        │   ├── metadata_merged_annotated.csv
        │   └── count/
        ├── B_VIC/
        │   ├── metadata_merged_annotated.csv
        │   └── count/
        ├── B_YAM/
        │   ├── metadata_merged_annotated.csv
        │   └── count/
        └── H5NX/
            ├── metadata_merged_annotated.csv
            └── count/
                ├── HA/
                ├── NA_N1/
                ├── NA_N2/
                └── ...
```

Count tables should be named with this pattern:

```text
aa_usage_by_<GROUPING>.csv
nt_usage_by_<GROUPING>.csv
```

Examples:

```text
aa_usage_by_HA_clade.csv
aa_usage_by_NA_clade.csv
aa_usage_by_HA_legacy_clade_yam.csv
aa_usage_by_Year_Month.csv
```

For H5NX-style datasets with multiple neuraminidase segments, keep each NA as its own gene folder, for example `count/NA_N1/`, `count/NA_N2/`, ..., `count/NA_N9/`. The app treats these as separate genes while preserving groupings such as `NA_subtype`, `HA_clade`, `Pathogenicity`, and `Year_Month`.

The app ignores validation-only count columns named `CodonStatus` and `CodonSource` when it builds caches. If a `Codon` column is present, it is preserved in the app cache as `Codon_Usage`.

### 5. Build or Refresh Caches

On first startup, FLUExplorer builds:

- `data/cache/FLU/app_cache_flu.rds`: compact metadata/statistics cache
- `data/cache/FLU/flu_explorer.duckdb`: normalized AA/NT usage table cache
- `data/cache/FLU/flu_explorer_duckdb_meta.rds`: DuckDB cache metadata

You can build/refresh them from R:

```r
Sys.setenv(FLUEXPLORER_REBUILD_FLU_CACHE = "true")
source("global.R", local = FALSE)
```

Unset `FLUEXPLORER_REBUILD_FLU_CACHE` after the rebuild if you do not want every startup to refresh the FLU caches.

```r
source("global.R")
```

The app checks raw CSV modification times and rebuilds stale caches when raw files under `data/raw` are newer than the existing cache files.

By default, DuckDB index creation is skipped to keep cache builds stable on large data. To opt in:

```r
Sys.setenv(FLUEXPLORER_DUCKDB_CREATE_INDEXES = "true")
source("global.R")
```

### 6. Run the App

From R:

```r
shiny::runApp(".")
```

Or from a shell in the repository:

```bash
Rscript -e 'shiny::runApp(".", host = "127.0.0.1", port = 4055)'
```

Then open:

```text
http://127.0.0.1:4055
```

## Repository Structure

```text
.
├── APP_INFO.md                  # Methods, reference datasets, and app update notes
├── DEVELOPMENT_LOG.md           # Development history and implementation notes
├── README.md                    # Setup and usage guide
├── global.R                     # Package loading, cache building, path helpers, query helpers
├── server.R                     # Shiny server logic and interactive analyses
├── ui.R                         # Shiny UI, navigation, and styling
├── data/
│   ├── raw/                     # User-provided raw metadata and count tables
│   ├── app_cache_flu.rds        # Generated compact metadata cache
│   ├── flu_explorer.duckdb      # Generated DuckDB usage cache
│   └── flu_explorer_duckdb_meta.rds
└── www/                         # Static app assets
```

## Troubleshooting

### The app starts slowly the first time

This is expected when caches are missing or stale. The first run reads raw CSVs and builds the compact metadata and DuckDB caches.

### A subtype or gene does not appear

Check that the subtype has:

- `data/raw/FLU/<SUBTYPE>/metadata_merged_annotated.csv`
- one or more count files under `data/raw/FLU/<SUBTYPE>/count/<GENE>/`
- count files named `aa_usage_by_*.csv` or `nt_usage_by_*.csv`

### NT mode has no genes

NT mode only appears when matching `nt_usage_by_*.csv` files exist. AA-only datasets will not expose NT gene choices.

### Memory use is high

Install `duckdb` and rebuild caches with `source("global.R")`. DuckDB-backed queries keep the app from loading full usage tables into R memory for most workflows.

## Authors

Lei Li - Initial work and development.

## License

This project is licensed under the MIT License.
