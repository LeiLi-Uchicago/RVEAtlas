# =============================================================================
# adapter-cache-build.R
#
# Build the "standard"-schema app cache (RSV, CHIKV, and any future adapter
# pathogen that uses schema = "standard") from raw output tables, so the app can
# regenerate a missing cache at startup instead of relying on an external
# pipeline. Mirrors the FLU builder's usage schema plus the Position_Key /
# Position_Label / Position_Offset columns the adapter needs.
#
# RAW INPUT (per pathogen under data/raw/<PATHOGEN>/):
#   Flat (single subtype, e.g. CHIKV):   output_tables/ + cleaned_metadata.tsv
#   Subtyped (e.g. RSV -> A, B):          <SUBTYPE>/output_tables/ + <SUBTYPE>/cleaned_metadata.tsv
#   output_tables/usage_<GENE>_by_clade.tsv       cols: position, aa, codon, count, clade, Year, Month
#   output_tables/usage_<GENE>_by_Year_month.tsv  cols: position, aa, codon, count, Year, Month
#   `position` looks like "F:262" (or "F:262+1" for an insertion).
#
# OUTPUT (matching PATHOGEN_ADAPTERS[[id]]$duckdb / $metadata):
#   data/cache/<PATHOGEN>/<name>_explorer.duckdb   -- `usage` table (standard schema)
#   data/cache/<PATHOGEN>/metadata_global.rds      -- per-sequence metadata for
#       Dataset Insights + Genetic Clade (conservation_entropy.rds and
#       dataset_insights.rds auto-generate at app start).
# =============================================================================

# Discover the build units for a standard pathogen: a list of
# list(group, out_dir, meta_tsv). `group` is the raw subtype value stored in the
# usage "Group" column (== adapter_subtype_value of the subtype choice).
sv_adapter_build_units <- function(pathogen_id, raw_root, cfg) {
  flat_out <- file.path(raw_root, "output_tables")
  if (dir.exists(flat_out)) {
    grp <- adapter_subtype_value(unname(cfg$subtype_choices)[[1]])
    return(list(list(group = grp, out_dir = flat_out,
                     meta_tsv = file.path(raw_root, "cleaned_metadata.tsv"))))
  }
  subdirs <- list.dirs(raw_root, recursive = FALSE, full.names = TRUE)
  units <- list()
  for (d in subdirs) {
    od <- file.path(d, "output_tables")
    if (dir.exists(od)) {
      units[[length(units) + 1]] <- list(
        group = basename(d), out_dir = od,
        meta_tsv = file.path(d, "cleaned_metadata.tsv"))
    }
  }
  units
}

# SQL that parses "F:262[+N]" into base/label/offset/key (Position_Key is the
# %09d__%03d form the adapter's position filter matches on).
.sv_pos_base   <- "CAST(split_part(split_part(position, ':', 2), '+', 1) AS INTEGER)"
.sv_pos_offset <- "CASE WHEN position LIKE '%+%' THEN CAST(split_part(position, '+', 2) AS INTEGER) ELSE 0 END"
.sv_pos_label  <- "split_part(position, ':', 2)"
.sv_pos_key    <- sprintf("printf('%%09d__%%03d', %s, %s)", .sv_pos_base, .sv_pos_offset)
# Map the raw deletion/missing markers onto the app's '-'/'X' convention.
.sv_aa_expr    <- "CASE WHEN aa = 'DEL' THEN '-' WHEN aa = 'MISSING' THEN 'X' ELSE aa END"
# Year-Month string, matching normalize_year_month_filter() (2-digit month).
.sv_ym_expr <- "
  CASE
    WHEN Year IS NULL OR trim(CAST(Year AS VARCHAR)) IN ('', 'NA') THEN NULL
    WHEN trim(CAST(Year AS VARCHAR)) IN ('Unknown','unassigned','Unassigned') THEN trim(CAST(Year AS VARCHAR))
    WHEN Month IS NULL OR trim(CAST(Month AS VARCHAR)) IN ('', 'NA', 'Unknown', 'unassigned', 'Unassigned')
      THEN trim(CAST(Year AS VARCHAR)) || '-Unknown'
    ELSE trim(CAST(Year AS VARCHAR)) || '-' ||
         CASE WHEN regexp_matches(trim(CAST(Month AS VARCHAR)), '^[0-9]+$')
              THEN lpad(trim(CAST(Month AS VARCHAR)), 2, '0')
              ELSE trim(CAST(Month AS VARCHAR)) END
  END"

.sv_read_tsv_sql <- function(path) {
  sprintf("read_csv('%s', delim='\t', header=true, all_varchar=true, ignore_errors=true, null_padding=true)",
          gsub("'", "''", path))
}

# gene from "usage_<GENE>_by_clade.tsv"; hyphens -> underscores (M2-1 -> M2_1).
.sv_gene_from_file <- function(fname, suffix) {
  g <- sub(paste0("^usage_(.+)_by_", suffix, "\\.tsv$"), "\\1", fname)
  gsub("-", "_", g)
}

sv_build_standard_usage_duckdb <- function(units, db_path) {
  tmp <- paste0(db_path, ".tmp")
  if (file.exists(tmp)) unlink(tmp)
  con <- DBI::dbConnect(duckdb::duckdb(), tmp, read_only = FALSE)
  on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  DBI::dbExecute(con, "PRAGMA memory_limit='3GB'")
  # threads=1: DuckDB's multi-threaded scan/index tear-down intermittently aborts
  # the R process on macOS ("mutex lock failed: Invalid argument"). Single-threaded
  # is a little slower but reliable for these one-off cache builds.
  DBI::dbExecute(con, "PRAGMA threads=1")
  DBI::dbExecute(con, '
    CREATE TABLE usage (
      "Group" VARCHAR, Variation_Type VARCHAR, Gene VARCHAR, Grouping_Type VARCHAR,
      Clade VARCHAR, Position INTEGER, Position_Label VARCHAR, Position_Key VARCHAR,
      Position_Offset INTEGER, AminoAcid VARCHAR, Count DOUBLE,
      Year VARCHAR, Month VARCHAR, Year_Month VARCHAR, Year_Month_Filter VARCHAR,
      Codon_Usage VARCHAR)')

  # Common projection of the position + aa + count columns (order matches table).
  proj <- sprintf(
    "%s AS Position, %s AS Position_Label, %s AS Position_Key, %s AS Position_Offset,
     %s AS AminoAcid, SUM(CAST(count AS DOUBLE)) AS Count",
    .sv_pos_base, .sv_pos_label, .sv_pos_key, .sv_pos_offset, .sv_aa_expr)

  for (u in units) {
    grp <- gsub("'", "''", u$group)
    clade_files <- list.files(u$out_dir, pattern = "^usage_.*_by_clade\\.tsv$")
    for (cf in clade_files) {
      gene <- .sv_gene_from_file(cf, "clade")
      cpath <- file.path(u$out_dir, cf)
      ympath <- file.path(u$out_dir, sub("_by_clade\\.tsv$", "_by_Year_month.tsv", cf))
      message("  ", u$group, " / ", gene)

      # Clade grouping (from by_clade)
      DBI::dbExecute(con, sprintf(
        "INSERT INTO usage
         SELECT '%s' AS \"Group\", 'AA' AS Variation_Type, '%s' AS Gene, 'Clade' AS Grouping_Type,
                CAST(clade AS VARCHAR) AS Clade, %s,
                NULL AS Year, NULL AS Month, NULL AS Year_Month, NULL AS Year_Month_Filter,
                ANY_VALUE(codon) AS Codon_Usage
         FROM %s GROUP BY ALL",
        grp, gene, proj, .sv_read_tsv_sql(cpath)))

      if (file.exists(ympath)) {
        # Year grouping (from by_Year_month, months collapsed)
        DBI::dbExecute(con, sprintf(
          "INSERT INTO usage
           SELECT '%s' AS \"Group\", 'AA' AS Variation_Type, '%s' AS Gene, 'Year' AS Grouping_Type,
                  CAST(Year AS VARCHAR) AS Clade, %s,
                  CAST(Year AS VARCHAR) AS Year, NULL AS Month, NULL AS Year_Month, NULL AS Year_Month_Filter,
                  ANY_VALUE(codon) AS Codon_Usage
           FROM %s GROUP BY ALL",
          grp, gene, proj, .sv_read_tsv_sql(ympath)))

        # Year_Month grouping (from by_Year_month)
        DBI::dbExecute(con, sprintf(
          "INSERT INTO usage
           SELECT '%s' AS \"Group\", 'AA' AS Variation_Type, '%s' AS Gene, 'Year_Month' AS Grouping_Type,
                  %s AS Clade, %s,
                  CAST(Year AS VARCHAR) AS Year, CAST(Month AS VARCHAR) AS Month,
                  %s AS Year_Month, %s AS Year_Month_Filter,
                  ANY_VALUE(codon) AS Codon_Usage
           FROM %s GROUP BY ALL",
          grp, gene, .sv_ym_expr, proj, .sv_ym_expr, .sv_ym_expr, .sv_read_tsv_sql(ympath)))
      }
    }
  }

  DBI::dbExecute(con, 'CREATE INDEX idx_usage_main ON usage ("Group", Gene, Grouping_Type, Position)')
  n <- DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM usage")$n
  DBI::dbDisconnect(con, shutdown = TRUE); on.exit()
  if (file.exists(db_path)) unlink(db_path)
  file.rename(tmp, db_path)
  message("  usage rows: ", format(n, big.mark = ","))
  invisible(n)
}

sv_build_standard_metadata <- function(units, meta_path) {
  # Read only the columns Dataset Insights / Genetic Clade need. A single
  # single-threaded connection is reused for every unit (opening/closing many
  # DuckDB connections can abort the R process on macOS -- see threads=1 note).
  want <- c("date", "region", "country", "host", "clade",
            "group_1", "group_2", "group_3", "group_4")
  con <- DBI::dbConnect(duckdb::duckdb(), ":memory:", read_only = FALSE)
  on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  DBI::dbExecute(con, "PRAGMA threads=1")
  frames <- lapply(units, function(u) {
    if (!file.exists(u$meta_tsv)) return(NULL)
    cols <- DBI::dbGetQuery(con, sprintf(
      "SELECT * FROM %s LIMIT 0", .sv_read_tsv_sql(u$meta_tsv)))
    present <- intersect(want, names(cols))
    sel <- paste(sprintf('CAST("%s" AS VARCHAR) AS "%s"', present, present), collapse = ", ")
    df <- DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s", sel, .sv_read_tsv_sql(u$meta_tsv)))
    df$Group <- u$group
    df
  })
  frames <- Filter(Negate(is.null), frames)
  meta <- dplyr::bind_rows(frames)
  if (nrow(meta) == 0) {
    warning("No metadata rows for ", meta_path)
  } else {
    # date is "YYYY-MM-DD" (or "XXXX-XX-XX" when unknown). Derive Year and a
    # year-month key named `YM` -- the column the Genetic Clade explorer reads
    # (metadata_global_to_summary derives YearMonth from it). Keep only valid
    # numeric year/month so unknown dates don't create junk month buckets.
    d  <- as.character(meta$date)
    yr <- substr(d, 1, 4)
    mo <- substr(d, 6, 7)
    valid_yr <- grepl("^[0-9]{4}$", yr)
    valid_mo <- grepl("^(0[1-9]|1[0-2])$", mo)
    meta$Year <- ifelse(valid_yr, yr, NA_character_)
    meta$YM   <- ifelse(valid_yr & valid_mo, paste0(yr, "-", mo), NA_character_)
  }
  n <- nrow(meta)
  dir.create(dirname(meta_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(
    metadata_global = meta,
    total_raw = format(n, big.mark = ","),
    total_parsed = format(n, big.mark = ","),
    built_at = Sys.time()
  ), meta_path)
  message("  metadata rows: ", format(n, big.mark = ","))
  invisible(n)
}

# COVID uses the "covid" schema; its builder already exists as a self-contained
# script. Run its duckdb + metadata phases (the "aux" phase re-sources global.R,
# so we deliberately skip it -- conservation/insights auto-generate at app start).
build_covid_adapter_cache <- function() {
  script <- file.path("tools", "build_covid_cache.R")
  if (!file.exists(script)) stop("COVID builder not found: ", script)
  for (phase in c("duckdb", "metadata")) {
    status <- system2("Rscript", c(shQuote(script), phase))
    if (!identical(as.integer(status), 0L)) {
      stop("build_covid_cache.R '", phase, "' phase failed (exit ", status, ")")
    }
  }
  invisible(TRUE)
}

# Build (or rebuild) the full standard-schema cache for one pathogen.
build_standard_adapter_cache <- function(pathogen_id) {
  cfg <- PATHOGEN_ADAPTERS[[pathogen_id]]
  if (is.null(cfg) || !identical(cfg$schema, "standard")) {
    stop("build_standard_adapter_cache: ", pathogen_id, " is not a standard-schema pathogen")
  }
  raw_root <- file.path("data", "raw", pathogen_id)
  units <- sv_adapter_build_units(pathogen_id, raw_root, cfg)
  if (length(units) == 0) stop("No raw output_tables found under ", raw_root)
  message("Building ", pathogen_id, " cache from raw (", length(units), " unit(s))...")
  dir.create(dirname(cfg$duckdb), recursive = TRUE, showWarnings = FALSE)
  sv_build_standard_usage_duckdb(units, cfg$duckdb)
  sv_build_standard_metadata(units, cfg$metadata)
  invisible(TRUE)
}
