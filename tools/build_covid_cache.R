#!/usr/bin/env Rscript
# =============================================================================
# build_covid_cache.R
#
# Rebuilds the COVID app cache from the raw output tables under data/raw/COVID/.
# Produces:
#   data/cache/COVID/covid_explorer.duckdb   -- the `usage` table (covid schema)
#   data/cache/COVID/metadata_summary.rds    -- Dataset Insights inputs
# (conservation_entropy.rds and dataset_insights.rds auto-generate at app start.)
#
# RAW INPUT (one pipeline / format shared across datasets):
#   data/raw/COVID/output_tables/usage_<Protein>_by_<GroupType>.tsv
#       cols: position, aa, codon, count, <GroupType>, Year, Month   (clade files)
#             position, aa, codon, count, Year, Month                 (Year_month file)
#       `position` looks like "S:614" or "S:614+1" (insertion).
#   data/raw/COVID/cleaned_metadata.tsv  -- per-sequence metadata (large).
#
# IMPORTANT -- Year-Month double-counting fix:
#   The upstream usage_*_by_Year_month.tsv files count every sequence TWICE
#   (their totals are exactly 2x the clade totals). We therefore IGNORE those
#   files and derive the Year_month grouping from a single full-coverage clade
#   file (Nextstrain_clade), which counts each sequence once. This is the same
#   principle the FLU builder uses (derive time groupings from one source).
#
#   Rscript tools/build_covid_cache.R
# =============================================================================
suppressWarnings(suppressMessages(library(DBI)))

RAW_DIR   <- "data/raw/COVID/output_tables"
META_TSV  <- "data/raw/COVID/cleaned_metadata.tsv"
OUT_DIR   <- "data/cache/COVID"
DUCKDB    <- file.path(OUT_DIR, "covid_explorer.duckdb")
META_RDS  <- file.path(OUT_DIR, "metadata_summary.rds")
SUBTYPE   <- "COVID:SARS-CoV-2"

# Clade grouping schemes present in the raw output tables. The first full-coverage
# one is also used to derive the Year_month grouping (see note above).
CLADE_SCHEMES <- c("Nextstrain_clade", "Nextclade_pango", "clade_who")
YM_SOURCE     <- "Nextstrain_clade"   # full coverage; used to derive Year_month

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

proteins_in_raw <- function() {
  files <- list.files(RAW_DIR, pattern = "^usage_.*_by_.*\\.tsv$")
  sort(unique(sub("^usage_(.+)_by_.*\\.tsv$", "\\1", files)))
}

raw_file <- function(protein, grouptype) {
  file.path(RAW_DIR, sprintf("usage_%s_by_%s.tsv", protein, grouptype))
}

# SQL expressions that parse the composite `position` ("S:614" / "S:614+1").
POS_EXPRS <- "
  position AS Position,
  split_part(position, ':', 2) AS PositionLabel,
  CAST(split_part(split_part(position, ':', 2), '+', 1) AS INTEGER) AS PositionBase,
  CASE WHEN position LIKE '%+%' THEN CAST(split_part(position, '+', 2) AS INTEGER) ELSE 0 END AS InsertionOffset,
  CAST(split_part(split_part(position, ':', 2), '+', 1) AS DOUBLE) AS PositionOrder"

# Normalize (Year, Month) -> "YYYY-MM" / "YYYY-Unknown", matching the app.
YM_EXPR <- "
  CASE
    WHEN Year IS NULL OR CAST(Year AS VARCHAR) IN ('', 'NA') THEN NULL
    WHEN CAST(Year AS VARCHAR) IN ('Unknown','unassigned','Unassigned') THEN CAST(Year AS VARCHAR)
    WHEN Month IS NULL OR CAST(Month AS VARCHAR) IN ('', 'NA', 'Unknown', 'unassigned', 'Unassigned')
      THEN CAST(Year AS VARCHAR) || '-Unknown'
    ELSE CAST(Year AS VARCHAR) || '-' || CAST(Month AS VARCHAR)
  END"

read_csv_sql <- function(path) {
  sprintf("read_csv('%s', delim='\t', header=true, all_varchar=true, ignore_errors=true)", path)
}

build_usage_duckdb <- function() {
  if (file.exists(DUCKDB)) file.remove(DUCKDB)
  con <- dbConnect(duckdb::duckdb(), DUCKDB, read_only = FALSE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, "PRAGMA memory_limit='4GB'")
  dbExecute(con, "PRAGMA threads=4")

  dbExecute(con, "CREATE TABLE usage (
    Protein VARCHAR, GroupType VARCHAR, GroupValue VARCHAR,
    Position VARCHAR, PositionLabel VARCHAR, PositionBase DOUBLE,
    InsertionOffset DOUBLE, PositionOrder DOUBLE,
    AA VARCHAR, Codon VARCHAR, Count DOUBLE)")

  proteins <- proteins_in_raw()
  message("Proteins: ", paste(proteins, collapse = ", "))

  for (p in proteins) {
    # --- clade groupings: one row-set per scheme, aggregated ---
    for (scheme in CLADE_SCHEMES) {
      f <- raw_file(p, scheme)
      if (!file.exists(f)) { message("  skip (missing): ", basename(f)); next }
      message("  ", p, " / ", scheme)
      dbExecute(con, sprintf(
        "INSERT INTO usage
         SELECT '%s' AS Protein, '%s' AS GroupType, CAST(\"%s\" AS VARCHAR) AS GroupValue,
                %s, aa AS AA, codon AS Codon, SUM(CAST(count AS DOUBLE)) AS Count
         FROM %s GROUP BY ALL",
        p, scheme, scheme, POS_EXPRS, read_csv_sql(f)))
    }
    # --- Year_month derived from the full-coverage clade file (no doubling) ---
    ym_src <- raw_file(p, YM_SOURCE)
    if (file.exists(ym_src)) {
      message("  ", p, " / Year_month (derived from ", YM_SOURCE, ")")
      dbExecute(con, sprintf(
        "INSERT INTO usage
         SELECT '%s' AS Protein, 'Year_month' AS GroupType, %s AS GroupValue,
                %s, aa AS AA, codon AS Codon, SUM(CAST(count AS DOUBLE)) AS Count
         FROM %s GROUP BY ALL",
        p, YM_EXPR, POS_EXPRS, read_csv_sql(ym_src)))
    }
  }
  n <- dbGetQuery(con, "SELECT COUNT(*) AS n FROM usage")$n
  message("usage rows: ", format(n, big.mark = ","))
}

build_metadata_summary <- function() {
  # Stage the (few) needed metadata columns into a temp DuckDB, then aggregate.
  stage <- tempfile(fileext = ".duckdb")
  con <- dbConnect(duckdb::duckdb(), stage, read_only = FALSE)
  on.exit({ dbDisconnect(con, shutdown = TRUE); unlink(stage) })
  dbExecute(con, "PRAGMA memory_limit='4GB'"); dbExecute(con, "PRAGMA threads=4")

  fill_cols <- c("clade_who", "Nextclade_pango", "Nextstrain_clade", "pango_lineage",
                 "region", "country", "division", "host")
  message("Staging metadata (this scans the large TSV once)...")
  dbExecute(con, sprintf(
    "CREATE TABLE meta AS SELECT
        CAST(Year AS VARCHAR) AS Year, %s AS YearMonth,
        %s
     FROM %s",
    YM_EXPR, paste(sprintf("\"%s\"", fill_cols), collapse = ", "), read_csv_sql(META_TSV)))
  total <- dbGetQuery(con, "SELECT COUNT(*) AS n FROM meta")$n
  message("metadata rows: ", format(total, big.mark = ","))

  # breakdowns: one aggregated df per (time_col, fill_col)
  breakdowns <- list()
  for (tcol in c("Year", "YearMonth")) {
    for (fcol in fill_cols) {
      df <- dbGetQuery(con, sprintf(
        "SELECT CAST(%s AS VARCHAR) AS XValue, CAST(\"%s\" AS VARCHAR) AS FillValue,
                COUNT(*) AS Count
         FROM meta WHERE %s IS NOT NULL AND CAST(%s AS VARCHAR) <> ''
         GROUP BY 1, 2", tcol, fcol, tcol, tcol))
      df$FillValue[is.na(df$FillValue) | df$FillValue == ""] <- "Unknown"
      df$XOrder <- df$XValue
      breakdowns[[paste(tcol, fcol, sep = "__")]] <- df
    }
  }

  time_plot <- dbGetQuery(con,
    "SELECT Year, COALESCE(NULLIF(clade_who,''),'Unknown') AS clade_who, COUNT(*) AS Count
     FROM meta WHERE Year IS NOT NULL AND Year <> '' GROUP BY 1, 2")
  region_plot <- dbGetQuery(con,
    "SELECT COALESCE(NULLIF(region,''),'Unknown') AS region, COUNT(*) AS Count FROM meta GROUP BY 1")
  gs <- dbGetQuery(con,
    "SELECT COUNT(*) AS total_sequences,
            COUNT(DISTINCT NULLIF(country,'')) AS countries_represented,
            MIN(TRY_CAST(Year AS INTEGER)) AS ymin, MAX(TRY_CAST(Year AS INTEGER)) AS ymax
     FROM meta")
  global_summary <- list(
    total_sequences = gs$total_sequences,
    countries_represented = gs$countries_represented,
    time_span = paste(gs$ymin, "-", gs$ymax))

  saveRDS(list(breakdowns = breakdowns, time_plot = time_plot,
               region_plot = region_plot, global_summary = global_summary), META_RDS)
  message("wrote ", META_RDS, " (", length(breakdowns), " breakdowns, ",
          format(global_summary$total_sequences, big.mark = ","), " sequences)")
}

args <- commandArgs(TRUE)
phase <- if (length(args) > 0) args[[1]] else "all"
if (phase %in% c("all", "duckdb"))   build_usage_duckdb()
if (phase %in% c("all", "metadata")) build_metadata_summary()
message("Done.")
