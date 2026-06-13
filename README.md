# RVEAtlas

![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)
![Shiny](https://img.shields.io/badge/Built_with-R_Shiny-success)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-purple)
![License](https://img.shields.io/badge/License-MIT-green)

**RVEAtlas** is an interactive Shiny application and companion project website for exploring respiratory virus evolution. The app supports clade-aware amino acid and nucleotide variation analysis across influenza, RSV, SARS-CoV-2, and a universal AAExplorer workflow for compatible NextAA outputs such as CHIKV.

The current UI includes a redesigned first page inspired by the `Learn/` reference design: a St. Jude-style teal navigation bar, blue/red/teal color system, a hero section, feature cards, full-width pathogen banner imagery, and a bottom research-context panel.

Project website: <https://leili-uchicago.github.io/RVEAtlas/>

AI setup skill: <https://github.com/LeiLi-Uchicago/RVEAtlas_Skill>

## Current Apps and Data Paths

- **FLU Explorer:** Human influenza amino acid and nucleotide variation across subtypes, clades, countries/regions, and collection dates.
- **RSV Explorer:** RSV A/B clade-aware amino acid variation.
- **COVID Explorer:** SARS-CoV-2 clade-aware amino acid variation.
- **AAExplorer:** Universal viewer for NextAA outputs from supported Nextclade3/Nextstrain datasets; current website example uses CHIKV.

The Shiny app uses the navbar pathogen and subtype controls to switch between available datasets. Dataset methods and notes are stored in `APP_INFO_FLU.md`, `APP_INFO_RSV.md`, `APP_INFO_COVID.md`, and `APP_INFO_CHIKV.md`.

## Key Features


### Dataset Insights

- Summarizes total sequences, represented countries, and collection time span.
- Visualizes sequencing volume over time, regional composition, and subtype or metadata breakdowns.

### Genetic Clade

- Lets users choose a clade annotation and search/select ranked clades.
- Shows sequence count, rank/share, active period, peak prevalence month, monthly prevalence curves, and metadata breakdowns.

### Single Position Explorer

- Shows amino acid or nucleotide distributions at a selected gene position.
- Supports grouping by year, year-month, clade, and available metadata annotations.
- Can display percentages or raw counts and export plots/tables.

### Pairwise Comparison

- Compares two selected groups or clades across genes.
- Identifies fixed or near-fixed differences using a user-defined dominant-frequency threshold.

### Gene-Wide Landscapes

- **Conservation:** Maps positional Shannon entropy across a gene.
- Variable-site shortcuts can jump directly into the Single Position Explorer.

## Website Updates

The static project website lives in `docs/` and is suitable for GitHub Pages.

Recent website updates include:

- A custom `docs/assets/favicon.svg` browser-tab icon.
- A highlighted **AI Setup** section for the RVEAtlas Skill.
- Explorer cards with icon+text actions for **Public user**, **St. Jude user**, and **Run locally**.
- Responsive typography and stat-card alignment improvements.

## AI-Assisted Local Setup

Users who want an AI coding agent to install and run one of the local Shiny apps can use:

<https://github.com/LeiLi-Uchicago/RVEAtlas_Skill>

Example prompt:

```text
Use the RVEAtlas skill in this repository to install and run RSVExplorer locally. Install it into ~/RVEAtlasApps. Start with a dry run, then do the real install if the dry run looks correct.
```

The skill includes deterministic installer logic for downloading app code, release data, missing R packages, and launching the selected Shiny app.

## Local Setup

### 1. Install R

Install R 4.0 or newer from [CRAN](https://cran.r-project.org/). RStudio is optional but useful for interactive development.

### 2. Clone the Repository

```bash
git clone https://github.com/LeiLi-Uchicago/RVEAtlas.git
cd RVEAtlas
```

### 3. Install R Packages

Open R from the project folder and install the required packages:

```r
install.packages(c(
  "shiny",
  "bslib",
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

`duckdb` and `DBI` are strongly recommended because the app uses compact DuckDB usage caches to avoid loading large count tables into memory. If `duckdb` is not installed, the app falls back to legacy RDS lazy loading when those cache files are available.

### 4. Organize Data

FLU raw metadata and count tables should be placed under `data/raw/FLU`. Each subtype should have its own folder:

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
        ├── B_VIC/
        ├── B_YAM/
        └── H5NX/
```

Count tables should be named with this pattern:

```text
aa_usage_by_<GROUPING>.csv
nt_usage_by_<GROUPING>.csv
```

For H5NX-style datasets with multiple neuraminidase segments, keep each NA as its own gene folder, for example `count/NA_N1/`, `count/NA_N2/`, ..., `count/NA_N9/`.

RSV, COVID, and CHIKV/AAExplorer data are expected as prebuilt cache assets under `data/cache/<PATHOGEN>/` according to the adapter paths in `global.R`.

### 5. Build or Refresh Caches

On first startup, the FLU workflow builds:

- `data/cache/FLU/app_cache_flu.rds`
- `data/cache/FLU/flu_explorer.duckdb`
- `data/cache/FLU/flu_explorer_duckdb_meta.rds`

You can force a rebuild from R:

```r
Sys.setenv(FLUEXPLORER_REBUILD_FLU_CACHE = "true")
source("global.R", local = FALSE)
```

Unset `FLUEXPLORER_REBUILD_FLU_CACHE` after the rebuild if you do not want every startup to refresh FLU caches.

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
├── APP_INFO_COVID.md            # COVID methods, source notes, and app information
├── APP_INFO_CHIKV.md            # AAExplorer/CHIKV methods and app information
├── APP_INFO_FLU.md              # FLU methods, source notes, and app information
├── APP_INFO_RSV.md              # RSV methods, source notes, and app information
├── DEVELOPMENT_LOG.md           # Development history and implementation notes
├── README.md                    # Setup and usage guide
├── global.R                     # Package loading, cache building, adapters, query helpers
├── server.R                     # Shiny server logic and interactive analyses
├── ui.R                         # Shiny UI, navigation, first-page design, and styling
├── docs/                        # Static GitHub Pages website
├── Learn/                       # Reference UI design used for the latest Home redesign
├── www/                         # Static app assets served by Shiny
└── data/
    ├── raw/                     # User-provided raw metadata and count tables
    └── cache/                   # Generated or prebuilt app caches
```

## Troubleshooting

### The app starts slowly the first time

This is expected when caches are missing or stale. The first run may read raw CSVs and build compact metadata and DuckDB caches.

### The navbar color or Home page does not update

Restart the Shiny app process, not only the browser tab. `bslib` and browser sessions can cache generated theme CSS during a running Shiny session.

### A subtype or gene does not appear

For FLU, check that the subtype has:

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
