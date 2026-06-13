# Development Log - FLU Divergence Explorer

## Date: May 11, 2026

### 1. Raw Data Layout Migration
- **New Raw Data Root:** Updated metadata and count-table discovery for the reorganized `data/raw/<subtype>/` layout.
- **Count Folder Support:** Updated count-table loading to read protein-specific files from `data/raw/<subtype>/count/<protein>/`.
- **Generated Cache Separation:** Added path helpers so generated legacy count `.rds` files, when needed, are written under `data/count_cache/` instead of back into the raw data folders.
- **Cache Freshness Checks:** Added raw CSV modification-time checks so compact RDS and DuckDB caches rebuild when reorganized or updated raw data are newer than existing caches.
- **Validation Column Exclusion:** Dropped `CodonStatus` and `CodonSource` before writing app caches while preserving `Codon` as the app-facing `Codon_Usage` field.
- **Mode-Aware Gene Discovery:** Restricted gene discovery by actual `aa_usage_by_*` or `nt_usage_by_*` files so AA-only raw data no longer exposes empty NT protein options.
- **Cache-Only Startup:** Updated startup and gene discovery so completed RDS/DuckDB caches are trusted when raw count files have been removed after cache generation.

---

## Date: May 5, 2026

### 1. Genetic Clade Explorer
- **New Top-Level Tab:** Added a Genetic Clade tab between Gene-Wide Landscapes and Single Position Explorer for subtype-specific clade metadata exploration.
- **Ranked Clade Selector:** Added clade annotation selection plus a searchable autocomplete clade/group selector with labels in the form `clade | rank #N | n=COUNT`.
- **Missing-Value Filtering:** Excluded blank, `Unknown`, unassigned, and other unknown-like metadata values from clade choices and ranked summaries.
- **Summary Cards:** Added cards for total sequences, rank/share within subtype, active period, and peak month/prevalence.
- **Monthly Prevalence View:** Added a Plotly line/area prevalence chart with light count bars on a secondary axis, monthly hover details, and a range slider for longer time series.
- **Metadata Breakdowns:** Added top country, region, and host breakdown plots plus a monthly detail table.

### 2. Compact Metadata Cache
- **Clade Summary Cache:** Extended `data/app_cache_flu.rds` with `metadata_clade_explorer`, storing compact summaries, monthly counts, subtype month totals, and top metadata breakdowns.
- **Automatic Refresh:** Older caches missing Genetic Clade summaries are rebuilt from raw `metadata_merged_annotated.csv` files when available.
- **Cache-Only Safety:** If a deployment has only an older cache and no raw metadata, the Genetic Clade tab shows a clear refresh notice instead of crashing.
- **Selectize Fix:** Updated the clade/group selector client-side for this input to preserve autocomplete/dropdown labels and avoid stale selected options when switching annotations.

---

## Date: April 28, 2026

### 1. DuckDB Data Backend
- **DuckDB Usage Cache:** Added a generated `data/flu_explorer.duckdb` cache that normalizes AA/NT usage CSVs into a single queryable `usage` table.
- **Query Pushdown:** Refactored Single Position, Pairwise Comparison, Entropy Landscape, and Mutation Lollipop workflows to request filtered/aggregated result sets from DuckDB rather than loading entire `.rds` tables into R memory.
- **Fallback Compatibility:** Preserved the legacy `.rds` lazy-loading path when the `duckdb` package or DuckDB cache is unavailable.
- **Index Guard:** Made DuckDB index creation opt-in via `FLUEXPLORER_DUCKDB_CREATE_INDEXES=true` after local testing showed index creation can be unstable on the full 421M-row cache.

### 2. Year-Month Filtering & Plot Fixes
- **Stable Time Filter UI:** Replaced the draggable Year-Month slider with explicit Start and End selectors in the Single Position Explorer to avoid jumpy behavior across dense time ranges.
- **Year-Month Normalization Fix:** Corrected DuckDB normalization so `Year_Month` grouping is derived from `Year` and `Month` instead of falling back to `Unknown`.
- **Backward-Compatible Query Fix:** Updated DuckDB query helpers to use `Year_Month_Filter` whenever `Grouping_Type == "Year_Month"`, so already-built caches with incorrect `Clade = Unknown` still return correct `YYYY-MM` groups.
- **Sparse Month Plot Fix:** Added safe Year-Month axis break sampling to prevent the `wrong sign in 'by' argument` error when a position has sparse or special time values.

---

## Date: March 26, 2026

### 1. Performance & User Experience (UX) Enhancements
- **Global & Modal Loaders:** Implemented full-screen `waiter` loading screens across the app (Global context switches, Single Position Explorer, and Pairwise Comparison) to provide immediate visual feedback and freeze the UI during heavy data processing.
- **Lazy Loading Tabs:** Refactored reactive observers to utilize `session$clientData` hidden states, ensuring that heavy computations only trigger when a specific tab is visible, eliminating background lag.
- **Removed Cross-Tab Syncing:** Decoupled input dependencies between tabs to prevent reactive cascades and double-loading glitches when changing groups or subtypes.

### 2. Memory Management (Posit Server Optimization)
- **Real-Time Memory Monitor:** Added a floating widget to the UI to track RAM usage in real-time, including a manual "Clear Cache" button to gracefully release memory.
- **Aggressive Garbage Collection:** Implemented automatic cache clearing and garbage collection (`gc()`) upon user session termination and during cache evictions.
- **Cache Tuning:** Reduced the `get_lazy_table()` LRU cache size from 5 to 3 tables and added startup memory flushes to keep the application footprint strictly under Posit Server's 1GB memory limit.

### 3. Feature Deprecation
- **MSA Tab Removal:** Commented out the Multiple Sequence Alignment (MSA) tab and its associated heavy Bioconductor dependencies (`msaR`, `Biostrings`, `msa`) to substantially conserve RAM and improve app initialization times.

---

## Date: March 16, 2026

### 1. Migration from RSV to FLU (H1N1/H3N2)
- **Multi-Subtype Data Pipeline:** Overhauled `global.R` to support dynamic loading of multiple influenza subtypes. The app now merges metadata and usage tables from `data/H1N1/` and `data/H3N2/` subdirectories.
- **Neuraminidase (NA) Protein Fix:** Implemented a critical fix in `read_csv` calls (`na = character()`) to prevent the "NA" protein name from being interpreted as a logical missing value.
- **Robust Column Mapping:** Standardized inconsistent naming conventions across source files (e.g., mapping `Protein` to `Gene` and `HA_clade` to `Clade`) using length-stable `rename_with` logic.

### 2. UI/UX Global Overhaul
- **Centralized Header Controls:** Moved the Subtype selector to a high-visibility, fixed-position container in the top-right header, alongside the Data Mode (AA/NT) switch.
- **Streamlined Workflow:** Removed redundant subtype selectors from all individual tabs, enabling a "select once, explore everywhere" workflow.
- **Branding Refresh:** Updated all interface labels, home page content, and documentation to reflect the focus on Influenza Virus diversity.
- **Safe Dropdown Handling:** Applied `na.omit()` and safe `1:nrow()` checks to all selectors to prevent Shiny crashes when handling empty or partially missing data frames (e.g., `important_pos_df`).

### 3. Data Integrity & Visualization
- **Subtype-Specific Coloring:** Implemented explicit color mapping for subtypes (H1N1: Blue, H3N2: Red) in sequencing stats and geographical plots.
- **Enhanced Cache Logic:** Incremented the RDS cache version (`v2`) to force a rebuild, ensuring the "NA" protein fix and merged multi-subtype structure are correctly applied.

---

## Date: March 6, 2026 (Legacy RSV Version)

### 1. Dual Variation Support (AA & NT)
- **Standardized Data Infrastructure:** Implemented a unified data loading pipeline in `global.R` that handles both Amino Acid (AA) and Nucleotide (NT) datasets.
- **Reactive Data Switcher:** Integrated a reactive backend in `server.R` that dynamically swaps data sources based on user selection.

... (rest of legacy logs)
