# RVEAtlas

![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)
![Shiny](https://img.shields.io/badge/Built_with-R_Shiny-success)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-purple)
![License](https://img.shields.io/badge/License-MIT-green)

**RVEAtlas** is an interactive Shiny application and companion project website for exploring respiratory virus evolution. The app supports clade-aware amino acid and nucleotide variation analysis across influenza, RSV, SARS-CoV-2, and a universal AAExplorer workflow for compatible NextAA outputs such as CHIKV.
    
Project website: <https://leili-uchicago.github.io/RVEAtlas/>

AI setup skill: <https://github.com/LeiLi-Uchicago/RVEAtlas_Skill>

## Current Datasets

- **FLU:** Human influenza amino acid and nucleotide variation across subtypes, clades, countries/regions, and collection dates.
- **RSV:** RSV A/B clade-aware amino acid variation.
- **COVID:** SARS-CoV-2 clade-aware amino acid variation.
- **CHIKV:** Universal viewer for NextAA outputs from supported Nextclade3/Nextstrain datasets; current website example uses CHIKV.

The Shiny app uses the navbar pathogen and subtype controls to switch between available datasets, with an animated loading overlay for smooth transitions. Only pathogens whose cache is present are selectable; the app can be deployed with any single pathogen (or any subset) and starts cleanly, with the others shown as "Unavailable". Dataset methods and notes are stored in `APP_INFO_FLU.md`, `APP_INFO_RSV.md`, `APP_INFO_COVID.md`, and `APP_INFO_CHIKV.md`. Shared platform-wide update notes are stored in `APP_INFO_PLATFORM.md` and displayed after each pathogen-specific Methods & Info page.

## Key Features

### Home / Pathogen Chooser

- Landing page with per-pathogen artwork cards, descriptions, and subtype chips, grouped into pathogen-specific and universal-configuration tiers.
- Cards mark the active dataset as "Viewing" and any pathogen without a cache as "Unavailable".

### Dataset Insights

- Summarizes total sequences, represented countries, and collection time span.
- Visualizes sequencing volume over time, regional composition, and subtype or metadata breakdowns.
- The Year / Year-Month breakdown supports a selectable sub-category (color) and highlights the COVID-19 pandemic window (Mar 2020 – May 2023) as a translucent band.

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
- **Structural/functional regions:** The entropy plot overlays translucent domain bands (e.g. HA1/HA2, Spike NTD/RBD/fusion peptide/HR1/HR2, RSV F2/p27/F1) defined in `www/structures/gene_regions.tsv`.
- Variable-site shortcuts can jump directly into the Single Position Explorer.

### 3D Protein Structure Views

Interactive `r3dmol` (3Dmol.js) structure views are available in Amino-Acid mode for the surface antigens of all four pathogens. Internal genes with no bundled structure hide the panel automatically.

- **Structures covered:**
  - **Influenza HA:** H1N1pdm09 (PDB 4JTV), seasonal H1N1 (6WCR), H3N2 (4HMG), and H5 (2FK0 clade 1, 9NRR clade 2.3.4.4b).
  - **Influenza NA:** N1 (3NSS) and N2 (4GZO) tetramers.
  - **RSV F:** postfusion trimer (3RKI), shared by RSV A and B.
  - **SARS-CoV-2 Spike:** closed prefusion trimer (6VXX).
- **Conservation on structure:** Colors residues by the currently-selected entropy basis (conserved → variable) with a gradient legend, so conserved cores and variable surfaces are visible in 3D.
- **Epitope + position mapping:** On the Single Site page, highlights curated antigenic sites for the active gene (e.g. H1 Sa/Sb/Ca1/Ca2/Cb, plus H3/H5/N1/N2, RSV-F, and Spike site sets) with a color-coded legend and marks the selected position in red with a residue label.
- **View styles:** Surface (default), Cartoon, Stick, or Sphere on both pages.
- **Config-driven:** Structures, chain maps, and epitope sets live in `www/structures/` (`structure_config.tsv` plus per-gene `epitopes_*.tsv`); bundled PDBs are renumbered to the app's HA/NA numbering (see `tools/renumber_h5.py`, `tools/renumber_na.py`). Additional pathogens/genes are added by dropping in a PDB and a config row.

### Year-Balanced Conservation

- The Conservation page defaults to a **year-balanced** Shannon entropy that weights each year equally, correcting for uneven sampling across years (which otherwise makes genuinely variable sites look conserved). A toggle switches back to the **all-sequence** basis.

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
  "duckdb",
  "r3dmol"
))
```

`duckdb` and `DBI` are strongly recommended because the app uses compact DuckDB usage caches to avoid loading large count tables into memory. If `duckdb` is not installed, the app falls back to legacy RDS lazy loading when those cache files are available.

### 4. Organize Data
Pre-built data can be downloaded from:
```text
https://github.com/LeiLi-Uchicago/RVEAtlas/releases/download/v0.9-alpha/CHIKV.zip
https://github.com/LeiLi-Uchicago/RVEAtlas/releases/download/v0.9-alpha/COVID.zip
https://github.com/LeiLi-Uchicago/RVEAtlas/releases/download/v0.9-alpha/FLU.zip
https://github.com/LeiLi-Uchicago/RVEAtlas/releases/download/v0.9-alpha/RSV.zip
```
Download datasets you're interested in, and unzip them under RVEAtlas/data/cache folder:
data/
└── cache/
    └── FLU/
    └── COVID/
    └── RSV/
    └── CHIKV/
Skip step 5 if you decided to use Pre-built data. 

### 5. Build Caches from raw data (optional)
FLU raw metadata and count tables should be placed under `data/raw/FLU`. Each subtype should have its own folder:

```text
data/
└── raw/
    └── FLU/
        ├── H1N1pdm09/
        │   ├── metadata_merged_annotated.csv
        │   └── count/
        │       ├── HA/
        │       │   ├── aa_usage_by_HA_clade.csv
        │       │   ├── aa_usage_by_NA_clade.csv
        │       │   └── aa_usage_by_Year_Month.csv
        │       └── ...
        ├── H1N1seasonal/
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

RSV, COVID, and CHIKV/AAExplorer are served through the adapter layer. Their caches under `data/cache/<PATHOGEN>/` (paths defined in `global.R`) can either be shipped prebuilt or built from raw NextAA-style outputs placed under `data/raw/<PATHOGEN>/` (see [Build or Refresh Caches](#5-build-or-refresh-caches)).

Each pathogen's cache is independent, so the app can be deployed with any one pathogen (e.g. only RSV) or any subset; pathogens without a cache are shown as "Unavailable" and cannot be selected.

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

For the adapter pathogens (RSV, COVID, CHIKV), any cache that is missing but has raw data present under `data/raw/<PATHOGEN>/` is built automatically at startup. You can also rebuild a pathogen's full cache (DuckDB + metadata + conservation/insights) offline in one step:

```bash
Rscript tools/build_pathogen_cache.R <FLU|RSV|COVID|CHIKV>
```

This handles all four pathogens (COVID is delegated to `tools/build_covid_cache.R`), so app launches stay fast instead of recomputing conservation entropy on first load. Conservation-entropy caches for every available pathogen are otherwise refreshed on startup as needed.

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
├── APP_INFO_PLATFORM.md         # Platform-wide update log shown for all pathogens
├── APP_INFO_RSV.md              # RSV methods, source notes, and app information
├── DEVELOPMENT_LOG.md           # Development history and implementation notes
├── README.md                    # Setup and usage guide
├── global.R                     # Package loading, cache building, adapters, query helpers
├── server.R                     # Shiny server logic and interactive analyses
├── ui.R                         # Shiny UI, navigation, first-page design, and styling
├── ha_numbering_review_table.csv # HA position → mature/H-numbering reference map
├── R/                           # Shared modules (conservation, structure view, adapter cache build)
├── tools/                       # Cache builders and PDB/region/epitope generators
├── docs/                        # Static GitHub Pages website
├── www/                         # Static app assets served by Shiny
│   └── structures/              # PDBs, structure_config.tsv, epitopes_*.tsv, gene_regions.tsv
└── data/
    ├── raw/                     # User-provided raw metadata and count tables
    └── cache/                   # Generated or prebuilt app caches (one folder per pathogen)
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
