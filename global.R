library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(readr)
library(tidyr)
library(openxlsx)
library(plotly)
# library(msaR) 
# library(Biostrings)
# library(msa)
library(waiter)
library(lubridate)
library(tidyverse)
# library(leaflet)             # Fixes: could not find function "leafletOutput"
# library(leaflet.minicharts)  # Required for the pie charts on the map
library(shinyWidgets)
library(shinyjs)

USE_DUCKDB <- requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)
if (!USE_DUCKDB) {
  message("Package 'duckdb' is not installed. Falling back to legacy RDS lazy loading.")
}

# Disable scientific notation for the session
options(scipen = 999)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x[[1]]))) y else x
}

# ==========================================
# 1. GLOBAL DATA LOADING & SETUP
# ==========================================
# Version 11: Read important_positions.csv
RDS_CACHE <- file.path("data", "cache", "FLU", "app_cache_flu.rds")
DUCKDB_CACHE <- file.path("data", "cache", "FLU", "flu_explorer.duckdb")
DUCKDB_META_CACHE <- file.path("data", "cache", "FLU", "flu_explorer_duckdb_meta.rds")
CACHE_SCHEMA_VERSION <- 5L
RAW_DATA_DIR <- file.path("data", "raw", "FLU")
COUNT_RDS_CACHE_DIR <- file.path("data", "cache", "FLU", "count_cache")
VALIDATION_ONLY_COUNT_COLS <- c("CodonStatus", "CodonSource")
FORCE_REBUILD_FLU_CACHE <- identical(tolower(Sys.getenv("FLUEXPLORER_REBUILD_FLU_CACHE", "false")), "true")

metadata_file_path <- function(subtype) {
  file.path(RAW_DATA_DIR, subtype, "metadata_merged_annotated.csv")
}

count_root_path <- function(subtype, var_type = "AA") {
  raw_count_root <- file.path(RAW_DATA_DIR, subtype, "count")
  legacy_root <- file.path("data", subtype, var_type)
  if (dir.exists(raw_count_root)) raw_count_root else legacy_root
}

count_cache_root_path <- function(subtype, var_type = "AA") {
  file.path(COUNT_RDS_CACHE_DIR, subtype, var_type)
}

count_cache_gene_path <- function(subtype, var_type, gene) {
  file.path(count_cache_root_path(subtype, var_type), gene)
}

count_cache_file_path <- function(subtype, var_type, gene, group_by) {
  file.path(count_cache_gene_path(subtype, var_type, gene), paste0(tolower(var_type), "_usage_by_", group_by, ".rds"))
}

raw_count_file_path <- function(subtype, var_type, gene, group_by) {
  file.path(count_root_path(subtype, var_type), gene, paste0(tolower(var_type), "_usage_by_", group_by, ".csv"))
}

sort_flu_subtypes <- function(subtypes) {
  subtypes <- unique(as.character(subtypes))
  c(sort(subtypes[grepl("^H", subtypes)]), sort(subtypes[!grepl("^H", subtypes)]))
}

discover_flu_subtypes <- function(raw_dir = RAW_DATA_DIR) {
  possible_dirs <- if (dir.exists(raw_dir)) {
    list.dirs(raw_dir, full.names = FALSE, recursive = FALSE)
  } else {
    character(0)
  }
  possible_dirs <- possible_dirs[!possible_dirs %in% c("", ".DS_Store")]
  has_metadata_file <- vapply(possible_dirs, function(d) file.exists(metadata_file_path(d)), logical(1))
  detected <- possible_dirs[has_metadata_file]
  if (length(detected) > 0) {
    sort_flu_subtypes(detected)
  } else {
    c("H1N1", "H3N2", "B_VIC", "B_YAM")
  }
}

count_gene_dirs <- function(subtype, var_type = "AA", prefer_cache = FALSE) {
  prefix <- paste0(tolower(var_type), "_usage_by_")
  roots <- if (isTRUE(prefer_cache)) {
    c(count_cache_root_path(subtype, var_type), count_root_path(subtype, var_type))
  } else {
    c(count_root_path(subtype, var_type), count_cache_root_path(subtype, var_type))
  }
  roots <- roots[dir.exists(roots)]
  if (length(roots) == 0) return(character(0))
  dirs <- unlist(lapply(roots, list.dirs, full.names = TRUE, recursive = FALSE), use.names = FALSE)
  dirs <- dirs[basename(dirs) != ".DS_Store"]
  dirs <- dirs[vapply(dirs, function(path) {
    length(list.files(path, pattern = paste0("^", prefix, ".*\\.(csv|rds)$"))) > 0
  }, logical(1))]
  dirs[!duplicated(basename(dirs))]
}

available_count_genes <- function(subtype, var_type = "AA", prefer_cache = FALSE) {
  sort(unique(basename(count_gene_dirs(subtype, var_type, prefer_cache = prefer_cache))))
}

latest_raw_data_mtime <- function() {
  files <- list.files(RAW_DATA_DIR, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(as.POSIXct(NA))
  max(file.info(files)$mtime, na.rm = TRUE)
}

flu_raw_data_signature <- function(subtypes = SUBTYPES) {
  raw_files <- list.files(RAW_DATA_DIR, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  raw_files <- raw_files[!grepl("/\\.DS_Store$", raw_files)]
  raw_files <- sort(raw_files)
  info <- file.info(raw_files)
  rel_files <- sub(paste0("^", normalizePath(RAW_DATA_DIR, mustWork = FALSE), "/?"), "", normalizePath(raw_files, mustWork = FALSE))

  list(
    subtypes = sort_flu_subtypes(subtypes),
    files = rel_files,
    sizes = as.numeric(info$size),
    mtimes = as.numeric(info$mtime)
  )
}

flu_raw_data_signature_matches <- function(signature, subtypes = SUBTYPES) {
  if (is.null(signature) || !is.list(signature)) return(FALSE)
  identical(signature, flu_raw_data_signature(subtypes))
}

cache_older_than_raw_data <- function(cache_path) {
  if (!file.exists(cache_path)) return(TRUE)
  raw_mtime <- latest_raw_data_mtime()
  if (is.na(raw_mtime)) return(FALSE)
  cache_mtime <- file.info(cache_path)$mtime
  is.na(cache_mtime) || cache_mtime < raw_mtime
}

strip_validation_count_cols <- function(df) {
  df[, setdiff(names(df), VALIDATION_ONLY_COUNT_COLS), drop = FALSE]
}

parse_metadata_year_month <- function(metadata) {
  year_values <- if ("Year" %in% names(metadata)) as.character(metadata$Year) else rep(NA_character_, nrow(metadata))
  month_values <- if ("Month" %in% names(metadata)) as.character(metadata$Month) else rep(NA_character_, nrow(metadata))

  parsed_year <- suppressWarnings(as.numeric(trimws(year_values)))
  month_clean <- trimws(month_values)
  month_clean[is.na(month_clean) | month_clean == ""] <- NA_character_
  month_clean <- ifelse(
    grepl("^[0-9]+$", month_clean),
    stringr::str_pad(month_clean, width = 2, side = "left", pad = "0"),
    month_clean
  )

  metadata$Year <- parsed_year
  metadata$Month <- month_clean
  metadata$YM <- normalize_year_month_filter(parsed_year, month_clean)
  metadata
}

empty_metadata_clade_explorer <- function() {
  list(
    month_totals = tibble::tibble(),
    summaries = tibble::tibble(),
    monthly = tibble::tibble(),
    breakdowns = tibble::tibble()
  )
}

is_unknown_metadata_value <- function(value) {
  text <- stringr::str_squish(as.character(value))
  is.na(text) | text == "" | stringr::str_to_lower(text) %in% c(
    "unknown", "na", "n/a", "none", "null", "?", "unassigned", "not assigned",
    "trace 0", "undetermined", "not determined"
  )
}

build_metadata_clade_explorer_summary <- function(metadata, annotation_cols) {
  annotation_cols <- annotation_cols[annotation_cols %in% names(metadata)]

  for (column in c(annotation_cols, "Group", "YM", "region", "country", "Host")) {
    if (!column %in% names(metadata)) metadata[[column]] <- NA_character_
  }

  month_totals <- metadata %>%
    filter(!is.na(.data$YM), .data$YM != "") %>%
    count(Group, YearMonth = .data$YM, name = "Total") %>%
    arrange(Group, YearMonth)

  if (length(annotation_cols) == 0 || nrow(metadata) == 0) {
    out <- empty_metadata_clade_explorer()
    out$month_totals <- month_totals
    return(out)
  }

  subtype_totals <- metadata %>%
    count(Group, name = "TotalSequences")

  base <- purrr::map_dfr(annotation_cols, function(annotation) {
    metadata %>%
      transmute(
        Group = as.character(.data$Group),
        Annotation = annotation,
        Clade = as.character(.data[[annotation]]),
        YearMonth = as.character(.data$YM),
        region = as.character(.data$region),
        country = as.character(.data$country),
        host = as.character(.data$Host)
      ) %>%
      filter(!is_unknown_metadata_value(.data$Clade))
  })

  if (nrow(base) == 0) {
    out <- empty_metadata_clade_explorer()
    out$month_totals <- month_totals
    return(out)
  }

  monthly <- base %>%
    filter(!is.na(.data$YearMonth), .data$YearMonth != "") %>%
    count(Group, Annotation, Clade, YearMonth, name = "Count") %>%
    left_join(month_totals, by = c("Group", "YearMonth")) %>%
    mutate(
      Total = dplyr::coalesce(.data$Total, 0L),
      Percent = dplyr::if_else(.data$Total > 0, (.data$Count / .data$Total) * 100, 0)
    ) %>%
    arrange(Group, Annotation, Clade, YearMonth)

  annotation_totals <- base %>%
    count(Group, Annotation, name = "AnnotatedTotal") %>%
    left_join(subtype_totals, by = "Group") %>%
    mutate(
      MissingAnnotationCount = pmax(.data$TotalSequences - .data$AnnotatedTotal, 0),
      AnnotationCoverage = dplyr::if_else(.data$TotalSequences > 0, (.data$AnnotatedTotal / .data$TotalSequences) * 100, 0)
    )

  totals <- base %>%
    count(Group, Annotation, Clade, name = "StrainCount") %>%
    left_join(annotation_totals, by = c("Group", "Annotation")) %>%
    arrange(Group, Annotation, desc(StrainCount), Clade) %>%
    group_by(Group, Annotation) %>%
    mutate(
      Rank = row_number(),
      DatasetShare = dplyr::if_else(.data$AnnotatedTotal > 0, (.data$StrainCount / .data$AnnotatedTotal) * 100, 0),
      TotalDatasetShare = dplyr::if_else(.data$TotalSequences > 0, (.data$StrainCount / .data$TotalSequences) * 100, 0)
    ) %>%
    ungroup()

  periods <- base %>%
    filter(!is.na(.data$YearMonth), .data$YearMonth != "") %>%
    group_by(Group, Annotation, Clade) %>%
    summarise(
      FirstMonth = min(.data$YearMonth),
      LastMonth = max(.data$YearMonth),
      .groups = "drop"
    )

  peaks <- monthly %>%
    group_by(Group, Annotation, Clade) %>%
    arrange(desc(.data$Percent), desc(.data$Count), .data$YearMonth, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      Group,
      Annotation,
      Clade,
      PeakMonth = YearMonth,
      PeakCount = Count,
      PeakPercent = Percent
    )

  summaries <- totals %>%
    left_join(periods, by = c("Group", "Annotation", "Clade")) %>%
    left_join(peaks, by = c("Group", "Annotation", "Clade")) %>%
    arrange(Group, Annotation, Rank)

  breakdowns <- purrr::map_dfr(c("country", "region", "host"), function(category) {
    base %>%
      filter(!is_unknown_metadata_value(.data[[category]])) %>%
      count(Group, Annotation, Clade, Category = category, Value = .data[[category]], name = "Count") %>%
      group_by(Group, Annotation, Clade, Category) %>%
      arrange(desc(.data$Count), .data$Value, .by_group = TRUE) %>%
      mutate(Rank = row_number()) %>%
      filter(.data$Rank <= 10) %>%
      ungroup()
  })

  list(
    month_totals = month_totals,
    summaries = summaries,
    monthly = monthly,
    breakdowns = breakdowns
  )
}

raw_metadata_available <- function(subtypes = SUBTYPES) {
  length(subtypes) > 0 && any(file.exists(vapply(subtypes, metadata_file_path, character(1))))
}

normalize_year_month_filter <- function(year, month) {
  year_chr <- trimws(as.character(year))
  month_chr <- trimws(as.character(month))
  month_chr <- ifelse(grepl("^[0-9]+$", month_chr), stringr::str_pad(month_chr, width = 2, side = "left", pad = "0"), month_chr)
  missing_month <- is.na(month_chr) | month_chr %in% c("", "NA", "Unknown", "unassigned", "Unassigned")

  dplyr::case_when(
    is.na(year_chr) | year_chr %in% c("", "NA") ~ NA_character_,
    year_chr %in% c("Unknown", "unassigned", "Unassigned") ~ year_chr,
    missing_month ~ paste0(year_chr, "-Unknown"),
    TRUE ~ paste0(year_chr, "-", month_chr)
  )
}

empty_usage_duckdb_df <- function() {
  data.frame(
    Group = character(),
    Variation_Type = character(),
    Gene = character(),
    Grouping_Type = character(),
    Clade = character(),
    Position = numeric(),
    AminoAcid = character(),
    Count = numeric(),
    Year = character(),
    Month = character(),
    Year_Month = character(),
    Year_Month_Filter = character(),
    Codon_Usage = character(),
    stringsAsFactors = FALSE
  )
}

normalize_usage_table <- function(df, subtype, var_type, gene_name, group_name) {
  df <- strip_validation_count_cols(df)
  df <- df %>%
    dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
    mutate(Group = subtype)

  if (var_type == "NT") {
    df <- df %>%
      dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide"))
  }

  if (!"Gene" %in% names(df)) df$Gene <- gene_name
  if (!"AminoAcid" %in% names(df)) df$AminoAcid <- NA_character_
  if (!"Count" %in% names(df)) df$Count <- NA_real_
  if (!"Position" %in% names(df)) df$Position <- NA_real_
  if (!"Year" %in% names(df)) df$Year <- NA_character_
  if (!"Month" %in% names(df)) df$Month <- NA_character_
  if (!"Year_Month" %in% names(df)) {
    if (group_name == "Year_Month") {
      df$Year_Month <- normalize_year_month_filter(df$Year, df$Month)
    } else {
      df$Year_Month <- NA_character_
    }
  }
  if (!group_name %in% names(df)) df[[group_name]] <- "Unknown"
  if (group_name == "Year_Month") df[[group_name]] <- normalize_year_month_filter(df$Year, df$Month)
  if (!"Codon_Usage" %in% names(df) && "Codon" %in% names(df)) df$Codon_Usage <- df$Codon
  if (!"Codon_Usage" %in% names(df)) df$Codon_Usage <- NA_character_

  data.frame(
    Group = as.character(df$Group),
    Variation_Type = as.character(var_type),
    Gene = as.character(df$Gene),
    Grouping_Type = as.character(group_name),
    Clade = as.character(df[[group_name]]),
    Position = suppressWarnings(as.numeric(df$Position)),
    AminoAcid = as.character(df$AminoAcid),
    Count = suppressWarnings(as.numeric(df$Count)),
    Year = as.character(df$Year),
    Month = as.character(df$Month),
    Year_Month = as.character(df$Year_Month),
    Year_Month_Filter = normalize_year_month_filter(df$Year, df$Month),
    Codon_Usage = as.character(df$Codon_Usage),
    stringsAsFactors = FALSE
  )
}

build_usage_duckdb_cache <- function(subtypes = SUBTYPES, db_path = DUCKDB_CACHE) {
  if (!USE_DUCKDB) return(FALSE)

  tmp_path <- paste0(db_path, ".tmp")
  if (file.exists(tmp_path)) unlink(tmp_path)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp_path, read_only = FALSE)
  on.exit({
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  DBI::dbExecute(con, "PRAGMA memory_limit='700MB'")
  DBI::dbWriteTable(con, "usage", empty_usage_duckdb_df(), overwrite = TRUE)

  usage_groups <- character()

  for (subtype in subtypes) {
    message("Processing usage tables into DuckDB for ", subtype)

    for (var_type in c("AA", "NT")) {
      var_root <- count_root_path(subtype, var_type)
      if (!dir.exists(var_root)) next

      gene_dirs <- count_gene_dirs(subtype, var_type)
      for (g_dir in gene_dirs) {
        gene_name <- basename(g_dir)
        message("  Processing ", var_type, " / ", gene_name)

        prefix <- if (var_type == "AA") "aa" else "nt"
        pattern <- paste0("^", prefix, "_usage_by_.*\\.csv$")
        files <- list.files(g_dir, pattern = pattern, full.names = TRUE)

        for (f in files) {
          group_name <- sub(paste0("^", prefix, "_usage_by_(.*)\\.csv$"), "\\1", basename(f))
          usage_groups <- unique(c(usage_groups, group_name))

          df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
            normalize_usage_table(subtype, var_type, gene_name, group_name)

          DBI::dbWriteTable(con, "usage", df, append = TRUE)

          if (group_name == "Year_Month") {
            year_df <- df %>%
              group_by(Group, Variation_Type, Gene, Position, AminoAcid, Year) %>%
              summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
              mutate(
                Grouping_Type = "Year",
                Clade = as.character(Year),
                Month = NA_character_,
                Year_Month = NA_character_,
                Year_Month_Filter = NA_character_,
                Codon_Usage = NA_character_
              ) %>%
              dplyr::select(names(empty_usage_duckdb_df()))

            DBI::dbWriteTable(con, "usage", year_df, append = TRUE)
            usage_groups <- unique(c(usage_groups, "Year"))
          }

          rm(df)
          gc(FALSE)
        }
      }
    }
  }

  if (identical(tolower(Sys.getenv("FLUEXPLORER_DUCKDB_CREATE_INDEXES", "false")), "true")) {
    DBI::dbExecute(con, "CREATE INDEX idx_usage_main ON usage (\"Group\", Variation_Type, Gene, Grouping_Type, Position)")
    DBI::dbExecute(con, "CREATE INDEX idx_usage_clade ON usage (\"Group\", Variation_Type, Gene, Grouping_Type, Clade)")
    DBI::dbExecute(con, "CREATE INDEX idx_usage_time ON usage (\"Group\", Variation_Type, Gene, Grouping_Type, Position, Year_Month_Filter)")
  } else {
    message("Skipping DuckDB index creation. Set FLUEXPLORER_DUCKDB_CREATE_INDEXES=true to enable it during cache builds.")
  }

  saveRDS(
    list(
      cache_schema_version = CACHE_SCHEMA_VERSION,
      usage_groups = usage_groups,
      built_at = Sys.time(),
      subtypes = sort_flu_subtypes(subtypes),
      raw_data_signature = flu_raw_data_signature(subtypes)
    ),
    DUCKDB_META_CACHE
  )
  DBI::dbDisconnect(con, shutdown = TRUE)

  if (file.exists(db_path)) unlink(db_path)
  file.rename(tmp_path, db_path)
}

ensure_usage_duckdb_cache <- function(force_rebuild = FORCE_REBUILD_FLU_CACHE) {
  if (!USE_DUCKDB) return(FALSE)
  if (file.exists(DUCKDB_CACHE) && !isTRUE(force_rebuild)) {
    if (file.exists(DUCKDB_META_CACHE)) {
      meta <- tryCatch(readRDS(DUCKDB_META_CACHE), error = function(e) NULL)
      meta_schema <- suppressWarnings(as.integer(meta$cache_schema_version %||% NA_integer_))
      if (!is.na(meta_schema) && identical(meta_schema, CACHE_SCHEMA_VERSION)) {
        return(TRUE)
      }
      message("DuckDB cache schema is outdated. Rebuilding from raw count tables...")
    } else {
      message("DuckDB cache metadata is missing. Rebuilding from raw count tables...")
    }
  }

  message("DuckDB cache missing or stale. Building: ", DUCKDB_CACHE)
  if (isTRUE(force_rebuild) && file.exists(DUCKDB_CACHE)) {
    message("FLUEXPLORER_REBUILD_FLU_CACHE=true. Rebuilding DuckDB usage cache from raw data.")
  }
  ok <- tryCatch(build_usage_duckdb_cache(), error = function(e) {
    message("DuckDB cache build failed: ", conditionMessage(e))
    FALSE
  })
  isTRUE(ok) && file.exists(DUCKDB_CACHE)
}

# Subtypes to load: dynamically detect valid subtype folders in data/raw/FLU.
# This includes H5NX and any future influenza subtype folder with metadata.
SUBTYPES <- discover_flu_subtypes()

if (file.exists(RDS_CACHE)) {
  # ---- FAST PATH: load everything from the pre-built cache ----
  message("Loading data from RDS cache: ", RDS_CACHE)
  cache <- readRDS(RDS_CACHE)
  
  # Check if the lightweight cache is present
  required_objects <- c("metadata_summary_stats", "important_pos_df", "metadata_grouping_cols")
  if (isTRUE(FORCE_REBUILD_FLU_CACHE) && raw_metadata_available()) {
    message("FLUEXPLORER_REBUILD_FLU_CACHE=true. Rebuilding FLU caches from raw data...")
    cache_loaded <- FALSE
  } else if (all(required_objects %in% names(cache))) {
    if (
      (is.null(cache$cache_schema_version) || !identical(as.integer(cache$cache_schema_version), CACHE_SCHEMA_VERSION)) &&
        raw_metadata_available()
    ) {
      message("RDS cache schema is outdated. Rebuilding from raw metadata...")
      cache_loaded <- FALSE
    } else if (!"metadata_clade_explorer" %in% names(cache) && raw_metadata_available()) {
      message("RDS cache is missing Genetic Clade summaries. Rebuilding from raw metadata...")
      cache_loaded <- FALSE
    } else {
      if (!is.null(cache$raw_data_signature) && !flu_raw_data_signature_matches(cache$raw_data_signature) && raw_metadata_available()) {
        message("RDS cache raw-data signature is stale. Set FLUEXPLORER_REBUILD_FLU_CACHE=true to rebuild from raw data.")
      } else if (is.null(cache$raw_data_signature) && raw_metadata_available()) {
        message("RDS cache has no raw-data signature. Set FLUEXPLORER_REBUILD_FLU_CACHE=true to rebuild and record the current raw layout.")
      }
      total_raw          <- cache$total_raw
      total_parsed       <- cache$total_parsed
      important_pos_df       <- cache$important_pos_df
      metadata_summary_stats <- cache$metadata_summary_stats
      total_countries_val    <- cache$total_countries_val
      time_range_val         <- cache$time_range_val
      metadata_groups        <- cache$metadata_groups
      metadata_groups <- sort_flu_subtypes(unique(c(metadata_groups, SUBTYPES)))
      metadata_years         <- cache$metadata_years
      metadata_grouping_cols <- cache$metadata_grouping_cols
      metadata_clade_explorer <- cache$metadata_clade_explorer
      if (is.null(metadata_clade_explorer)) {
        metadata_clade_explorer <- empty_metadata_clade_explorer()
      }

      rm(cache)   # free the wrapper list from memory
      message("RDS cache loaded successfully.")
      cache_loaded <- TRUE
    }
  } else {
    message("RDS cache is outdated. Rebuilding...")
    cache_loaded <- FALSE
  }
} else {
  cache_loaded <- FALSE
}

if (!cache_loaded) {
  # ---- SLOW PATH: process raw CSVs and build cache ----
  message("RDS cache not found. Processing raw CSV files...")
  
  # --- Metadata ---
  all_metadata <- list()
  for (subtype in SUBTYPES) {
    meta_path <- metadata_file_path(subtype)
    if (file.exists(meta_path)) {
      message("Loading metadata for ", subtype)
      # CRITICAL: na = character() ensures 'NA' (Neuraminidase) is NOT treated as a missing value
      meta <- read_csv(meta_path, col_types = cols(.default = "c"), show_col_types = FALSE, na = character())
      
      # Rename columns to match app expectations
      # Group = Subtype (A / H1N1)
      meta <- meta %>% 
        dplyr::rename(any_of(c(Group = "Subtype",
                               clade = "HA_clade",
                               G_clade = "NA_clade")))
      
      # Force Group name to match the folder exactly (bypassing messy fasta labels like "B / Victoria")
      meta$Group <- subtype
      
      # Parse Location into region and country
      loc_split <- stringr::str_split(meta$Location, " / ", simplify = TRUE)
      meta$region <- loc_split[, 1]
      meta$country <- loc_split[, 2]
      
      all_metadata[[subtype]] <- meta
    }
  }
  metadata_global <- bind_rows(all_metadata)
  
  # Identify dynamic grouping columns by matching with actual usage file groupings
  usage_groups <- c()
  for (subtype in SUBTYPES) {
    aa_root <- count_root_path(subtype, "AA")
    if (dir.exists(aa_root)) {
      # Look into the first gene subdirectory found
      gene_dirs <- count_gene_dirs(subtype, "AA")
      if (length(gene_dirs) > 0) {
        files <- list.files(gene_dirs[1], pattern = "^aa_usage_by_.*\\.csv")
        usage_groups <- unique(c(usage_groups, sub("^aa_usage_by_(.*)\\.csv$", "\\1", files)))
      }
    }
  }
  # Add Year explicitly since we will generate it from Year_Month
  if ("Year_Month" %in% usage_groups) {
    usage_groups <- unique(c(usage_groups, "Year"))
  }
  mapped_usage_groups <- usage_groups
  mapped_usage_groups[mapped_usage_groups == "HA_clade"] <- "clade"
  mapped_usage_groups[mapped_usage_groups == "NA_clade"] <- "G_clade"
  
  metadata_grouping_cols <- intersect(colnames(metadata_global), setdiff(mapped_usage_groups, c("Year", "Year_Month")))
  
  # Clean up empty clade names
  for(col in metadata_grouping_cols) {
    if(col %in% colnames(metadata_global)) {
      metadata_global[[col]] <- ifelse(is.na(metadata_global[[col]]) | metadata_global[[col]] == "" | metadata_global[[col]] == "trace 0", "Unknown", metadata_global[[col]])
    }
  }

  # Metadata time handling: use the explicit parsed Year/Month columns from
  # metadata_merged_annotated.csv. Do not parse Collection_Date here.
  metadata_global <- parse_metadata_year_month(metadata_global)
  
  total_raw    <- scales::comma(nrow(metadata_global))
  metadata_global <- metadata_global %>% filter(!is.na(Year))
  total_parsed <- scales::comma(nrow(metadata_global))
  
  # --- Process Usage Tables ---
  duckdb_ready <- ensure_usage_duckdb_cache()
  if (!duckdb_ready) {
    message("DuckDB is unavailable. Processing usage tables into legacy RDS files.")

    for (subtype in SUBTYPES) {
      message("Processing usage tables into RDS for ", subtype)

      for (var_type in c("AA", "NT")) {
        var_root <- count_root_path(subtype, var_type)
        if (!dir.exists(var_root)) next

        gene_dirs <- count_gene_dirs(subtype, var_type)
        for (g_dir in gene_dirs) {
          gene_name <- basename(g_dir)
          message("  Processing ", var_type, " / ", gene_name)

          prefix <- if(var_type == "AA") "aa" else "nt"
          pattern <- paste0("^", prefix, "_usage_by_.*\\.csv$")
          files <- list.files(g_dir, pattern = pattern, full.names = TRUE)

          for (f in files) {
            is_ym <- grepl(paste0(prefix, "_usage_by_Year_Month\\.csv$"), f)
            group_name <- sub(paste0("^", prefix, "_usage_by_(.*)\\.csv$"), "\\1", basename(f))

            df <- read_csv(f, show_col_types = FALSE, na = character()) %>%
              strip_validation_count_cols() %>%
              dplyr::rename_with(~ gsub("^Protein$", "Gene", .x), any_of("Protein")) %>%
              mutate(Group = subtype)

            if (var_type == "NT") {
              df <- df %>%
                dplyr::rename_with(~ gsub("^Nucleotide$", "AminoAcid", .x), any_of("Nucleotide"))
            }

            if (is_ym) {
              df <- df %>%
                mutate(Year_Month = normalize_year_month_filter(Year, Month))
            }

            if (!"Codon_Usage" %in% names(df) && "Codon" %in% names(df)) df$Codon_Usage <- df$Codon

            out_dir <- count_cache_gene_path(subtype, var_type, gene_name)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            saveRDS(df, file.path(out_dir, paste0(prefix, "_usage_by_", group_name, ".rds")))

            if (is_ym) {
              year_df <- df %>%
                group_by(Group, Gene, Position, Year, AminoAcid) %>%
                summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

              saveRDS(year_df, file.path(out_dir, paste0(prefix, "_usage_by_Year.rds")))
            }
          }
        }
      }
    }
  }

  # --- Important positions ---
  important_pos_file <- "data/important_positions.csv"
  if (file.exists(important_pos_file)) {
    message("Loading important positions from CSV...")
    important_pos_df <- read_csv(important_pos_file, show_col_types = FALSE)
    if (!"label" %in% names(important_pos_df)) {
      important_pos_df$label <- paste(important_pos_df$Subtype, important_pos_df$Gene, important_pos_df$Position, sep = " - ")
    }
  } else {
    important_pos_df <- data.frame(
      Gene = character(),
      Subtype = character(),
      Position = numeric(),
      Mutation = character(),
      Epitope = character(),
      Clinical_impact = character(),
      Source = character(),
      label = character(),
      stringsAsFactors = FALSE
    )
  }
  
  # --- Save RDS cache ---
  message("Writing RDS cache to: ", RDS_CACHE)
  
  # Pre-calculate these before saving if they weren't loaded from cache
  message("Pre-aggregating metadata statistics...")
  metadata_summary_stats <- metadata_global %>%
    mutate(YearMonth = .data$YM) %>%
    group_by(across(all_of(c("Year", "YearMonth", "Group", "region", "country", metadata_grouping_cols)))) %>%
    summarise(n = n(), .groups = "drop")

  metadata_clade_explorer <- build_metadata_clade_explorer_summary(metadata_global, metadata_grouping_cols)
  
  total_countries_val <- length(unique(metadata_global$country))
  time_range_val <- paste(min(metadata_global$Year, na.rm=T), "-", max(metadata_global$Year, na.rm=T))
  
  raw_groups <- na.omit(unique(metadata_global$Group))
  metadata_groups <- c(sort(raw_groups[grepl("^H", raw_groups)]), sort(raw_groups[!grepl("^H", raw_groups)]))
  
  metadata_years <- sort(na.omit(unique(metadata_global$Year)), decreasing = TRUE)

  saveRDS(
    list(
      cache_schema_version   = CACHE_SCHEMA_VERSION,
      raw_data_signature     = flu_raw_data_signature(SUBTYPES),
      total_raw              = total_raw,
      total_parsed           = total_parsed,
      important_pos_df       = important_pos_df,
      metadata_summary_stats = metadata_summary_stats,
      total_countries_val    = total_countries_val,
      time_range_val         = time_range_val,
      metadata_groups        = metadata_groups,
      metadata_years         = metadata_years,
      metadata_grouping_cols = metadata_grouping_cols,
      metadata_clade_explorer = metadata_clade_explorer
    ),
    file = RDS_CACHE
  )
  
  # --- STARTUP MEMORY FLUSH ---
  suppressWarnings(rm(all_metadata, metadata_global, meta))
  gc(verbose = FALSE)
}

if (!exists("metadata_clade_explorer") || is.null(metadata_clade_explorer)) {
  metadata_clade_explorer <- empty_metadata_clade_explorer()
}

duckdb_cache_ready <- ensure_usage_duckdb_cache(force_rebuild = FALSE)

# ==========================================
# 2. COORDINATE LOOKUP DATA (Pre-calculated for Performance)
# ==========================================
# region_coords <- data.frame(
#   region = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania"),
#   lat = c(1.0, 34.0, 48.0, 45.0, -15.0, -25.0),
#   lng = c(17.0, 100.0, 10.0, -100.0, -60.0, 135.0),
#   stringsAsFactors = FALSE
# )
# 
# # Move world_coords out of server.R to global.R to avoid recalculation on every session
# # ggplot2::map_data is slow, so we do it once here.
# message("Pre-calculating world coordinates...")
# world_coords <- ggplot2::map_data("world") %>%
#   dplyr::group_by(region) %>%
#   dplyr::summarise(lat = mean(lat), lng = mean(long), .groups = "drop") %>%
#   dplyr::rename(country = region)

# ---- Post-load steps ----
ALL_AAS <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*","X", "-")

aa_colors <- c(
  "A"="#E41A1C", "C"="#377EB8", "D"="#4DAF4A", "E"="#984EA3",
  "F"="#FF7F00", "G"="#FFFF33", "H"="#A65628", "I"="#F781BF",
  "K"="#999999", "L"="#66C2A5", "M"="#FC8D62", "N"="#8DA0CB",
  "P"="#E78AC3", "Q"="#A6D854", "R"="#FFD92F", "S"="#E5C494",
  "T"="#B3B3B3", "V"="#1B9E77", "W"="#D95F02", "Y"="#7570B3",
  "*"="#000000", "X"="#D3D3D3", "-"="#808080"
)

nt_colors <- c(
  "a"="#E41A1C", "c"="#377EB8", "g"="#4DAF4A", "t"="#984EA3",
  "A"="#E41A1C", "C"="#377EB8", "G"="#4DAF4A", "T"="#984EA3",
  "N"="#000000", "n"="#000000", "-"="#808080"
)

ggmsa_custom_colors <- data.frame(
  names = names(aa_colors),
  color = unname(aa_colors),
  stringsAsFactors = FALSE
)

# Generate rainbow palettes for all possible clades found in any clade column
if (!exists("metadata_grouping_cols")) {
  metadata_grouping_cols <- c("clade", "G_clade") # fallback
}

all_possible_clades <- unique(unlist(lapply(metadata_grouping_cols, function(col) metadata_summary_stats[[col]])))
all_possible_clades <- sort(na.omit(all_possible_clades))

# Exclude 'Unknown' from the main rainbow generation to give it a neutral color
actual_clades <- setdiff(all_possible_clades, "Unknown")
master_clade_colors <- viridis::viridis(length(actual_clades))
names(master_clade_colors) <- actual_clades

# Add neutral color for Unknown
master_clade_colors["Unknown"] <- "#d3d3d3"

# For backward compatibility with server logic that expects specific names
clade_colors_vec   <- master_clade_colors
g_clade_colors_vec <- master_clade_colors

# ==========================================
# 3. LAZY-LOADING CACHE MECHANISM (LRU)
# ==========================================
lazy_cache <- new.env(parent = emptyenv())
lazy_cache$keys <- character(0)
lazy_cache$data <- list()

LAZY_CACHE_MAX_TABLES <- 2
LAZY_CACHE_MAX_MEM_MB <- 450

get_lazy_table <- function(rds_path, max_tables = LAZY_CACHE_MAX_TABLES, max_mem_mb = LAZY_CACHE_MAX_MEM_MB) {
  if (!file.exists(rds_path)) return(NULL)
  
  if (rds_path %in% lazy_cache$keys) {
    # Move to the end of the line (most recently used)
    lazy_cache$keys <- c(setdiff(lazy_cache$keys, rds_path), rds_path)
    return(lazy_cache$data[[rds_path]])
  }
  
  # Check memory usage against threshold before attempting to load new data
  if (sum(gc(verbose = FALSE)[, 2]) > max_mem_mb) {
    message("Memory usage exceeds ", max_mem_mb, " MB threshold. Clearing cache...")
    lazy_cache$keys <- character(0)
    lazy_cache$data <- list()
    gc(verbose = FALSE)
  }

  # Read from disk and store in cache
  df <- readRDS(rds_path)
  lazy_cache$keys <- c(lazy_cache$keys, rds_path)
  lazy_cache$data[[rds_path]] <- df
  
  # Evict oldest if limit exceeded
  while (length(lazy_cache$keys) > max_tables) {
    evict <- lazy_cache$keys[1]
    lazy_cache$keys <- lazy_cache$keys[-1]
    lazy_cache$data[[evict]] <- NULL
    
    # Force memory release back to the OS
    gc(verbose = FALSE)
  }
  
  return(df)
}

# ==========================================
# 4. DUCKDB QUERY HELPERS
# ==========================================
usage_db_env <- new.env(parent = emptyenv())
usage_db_env$con <- NULL

usage_duckdb_available <- function() {
  isTRUE(USE_DUCKDB) && isTRUE(duckdb_cache_ready) && file.exists(DUCKDB_CACHE)
}

usage_db_conn <- function() {
  if (!usage_duckdb_available()) return(NULL)

  if (is.null(usage_db_env$con) || !DBI::dbIsValid(usage_db_env$con)) {
    usage_db_env$con <- DBI::dbConnect(duckdb::duckdb(), dbdir = DUCKDB_CACHE, read_only = TRUE)
    DBI::dbExecute(usage_db_env$con, "PRAGMA memory_limit='700MB'")
  }

  usage_db_env$con
}

usage_query <- function(sql, params = NULL) {
  con <- usage_db_conn()
  if (is.null(con)) return(NULL)
  DBI::dbGetQuery(con, sql, params = params)
}

usage_sql_in_values <- function(values) {
  con <- usage_db_conn()
  if (is.null(con) || length(values) == 0) return(NULL)
  paste(DBI::dbQuoteString(con, values), collapse = ", ")
}

usage_file_groups <- function(subtype, var_type, gene) {
  dirs <- c(
    count_cache_gene_path(subtype, var_type, gene),
    file.path(count_root_path(subtype, var_type), gene)
  )
  files <- unlist(lapply(dirs[dir.exists(dirs)], function(dir_path) {
    list.files(dir_path, pattern = paste0("^", tolower(var_type), "_usage_by_.*\\.(rds|csv)$"))
  }), use.names = FALSE)
  sort(unique(sub(paste0("^", tolower(var_type), "_usage_by_(.*)\\.(rds|csv)$"), "\\1", files)))
}

usage_available_genes <- function(subtype, var_type) {
  if (!usage_duckdb_available()) return(available_count_genes(subtype, var_type, prefer_cache = TRUE))

  res <- usage_query(
    "SELECT DISTINCT Gene FROM usage WHERE \"Group\" = ? AND Variation_Type = ?",
    list(subtype, var_type)
  )
  if (is.null(res) || nrow(res) == 0) return(available_count_genes(subtype, var_type, prefer_cache = TRUE))
  sort(stats::na.omit(as.character(res$Gene)))
}

usage_available_groups <- function(subtype, var_type, gene) {
  if (!usage_duckdb_available()) return(usage_file_groups(subtype, var_type, gene))

  res <- usage_query(
    "SELECT DISTINCT Grouping_Type FROM usage WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ?",
    list(subtype, var_type, gene)
  )
  groups <- sort(as.character(res$Grouping_Type))

  if (length(groups) == 0) usage_file_groups(subtype, var_type, gene) else groups
}

usage_distinct_group_values <- function(subtype, var_type, gene, group_by) {
  res <- usage_query(
    "SELECT DISTINCT
       CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade
     FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?",
    list(subtype, var_type, gene, group_by)
  )
  if (is.null(res) || nrow(res) == 0) return(character(0))

  values <- sort(stats::na.omit(as.character(res$Clade)))
  special_values <- c("Unknown", "unassigned", "Unassigned")
  present_specials <- intersect(special_values, values)
  if (length(present_specials) > 0) values <- c(setdiff(values, present_specials), present_specials)
  values
}

usage_max_position <- function(subtype, var_type, gene) {
  res <- usage_query(
    "SELECT MAX(Position) AS max_position FROM usage WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ?",
    list(subtype, var_type, gene)
  )
  if (is.null(res) || nrow(res) == 0 || is.na(res$max_position[1])) return(NA_real_)
  as.numeric(res$max_position[1])
}

usage_year_month_choices <- function(subtype, var_type, gene, group_by, position) {
  res <- usage_query(
    "SELECT DISTINCT Year_Month_Filter FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
       AND Position = ? AND Year_Month_Filter IS NOT NULL",
    list(subtype, var_type, gene, group_by, as.numeric(position))
  )
  if (is.null(res) || nrow(res) == 0) return(character(0))

  ym_values <- stats::na.omit(as.character(res$Year_Month_Filter))
  special_values <- c("Unknown", "unassigned", "Unassigned")
  present_specials <- intersect(special_values, ym_values)
  chronological_yms <- sort(setdiff(ym_values, special_values))
  c(present_specials, chronological_yms)
}

flu_position_choices <- function(subtype, var_type, gene) {
  if (usage_duckdb_available()) {
    res <- usage_query(
      "SELECT DISTINCT Position FROM usage WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? ORDER BY Position",
      list(subtype, var_type, gene)
    )
    if (!is.null(res) && nrow(res) > 0) {
      labels <- as.character(res$Position)
      return(stats::setNames(labels, labels))
    }
  }
  character(0)
}

usage_single_position <- function(subtype, var_type, gene, group_by, position, allowed_yms = NULL, min_seqs = 1, hide_empty_years = FALSE) {
  ym_filter <- ""
  if (!is.null(allowed_yms) && length(allowed_yms) > 0) {
    in_values <- usage_sql_in_values(allowed_yms)
    if (!is.null(in_values)) {
      ym_filter <- paste0(" AND Year_Month_Filter IN (", in_values, ")")
    }
  }

  sql <- paste0(
    "SELECT \"Group\", Gene, Position,
            CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
            AminoAcid, SUM(Count) AS Count,
            ANY_VALUE(Codon_Usage) AS Codon_Usage
     FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
       AND Position = ? AND AminoAcid NOT IN ('X', '-')",
    ym_filter,
    " GROUP BY \"Group\", Gene, Position,
       CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
       AminoAcid"
  )
  res <- usage_query(sql, list(subtype, var_type, gene, group_by, as.numeric(position)))
  if (is.null(res)) return(NULL)
  if (nrow(res) == 0) {
    out <- data.frame(Group=character(), Gene=character(), Position=numeric(), AminoAcid=character(), Count=numeric(), Valid_Total=numeric(), `Frequency(%)`=numeric(), check.names = FALSE)
    out[[group_by]] <- character()
    return(out)
  }

  res <- res %>%
    dplyr::rename(!!group_by := Clade) %>%
    group_by(.data[[group_by]]) %>%
    mutate(
      Valid_Total = sum(Count, na.rm = TRUE),
      `Frequency(%)` = (Count / Valid_Total) * 100
    ) %>%
    ungroup() %>%
    filter(Valid_Total >= min_seqs)

  if (group_by == "Year" && isTRUE(hide_empty_years)) {
    res <- res %>% filter(Valid_Total > 0)
  }

  if (all(is.na(res$Codon_Usage))) {
    res$Codon_Usage <- NULL
  }

  res
}

usage_pairwise_gene_data <- function(subtype, var_type, gene, group_by, clades = NULL) {
  clade_filter <- ""
  if (!is.null(clades) && length(clades) > 0) {
    in_values <- usage_sql_in_values(clades)
    if (!is.null(in_values)) {
      clade_filter <- paste0(
        " AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END IN (",
        in_values,
        ")"
      )
    }
  }

  sql <- paste0(
    "SELECT \"Group\", Gene,
            CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
            Position, AminoAcid, SUM(Count) AS Count,
            ANY_VALUE(Codon_Usage) AS Codon_Usage
     FROM usage
     WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?",
    clade_filter,
    " GROUP BY \"Group\", Gene,
       CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
       Position, AminoAcid"
  )
  res <- usage_query(sql, list(subtype, var_type, gene, group_by))
  if (is.null(res)) return(NULL)
  if (all(is.na(res$Codon_Usage))) res$Codon_Usage <- NULL
  res
}

usage_pairwise_differences_for_gene <- function(subtype, var_type, gene, group_by, clade1, clade2, min_freq) {
  res <- usage_query(
    "WITH agg AS (
       SELECT Gene, Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
         AminoAcid, SUM(Count) AS Variant_Count
       FROM usage
       WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
         AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END IN (?, ?)
         AND AminoAcid NOT IN ('X', '-')
       GROUP BY Gene, Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
         AminoAcid
     ),
     freq AS (
       SELECT *,
         SUM(Variant_Count) OVER (PARTITION BY Gene, Position, Clade) AS Total_Seqs,
         100.0 * Variant_Count / SUM(Variant_Count) OVER (PARTITION BY Gene, Position, Clade) AS Freq
       FROM agg
     ),
     ranked AS (
       SELECT *,
         ROW_NUMBER() OVER (PARTITION BY Gene, Position, Clade ORDER BY Freq DESC, AminoAcid) AS rn
       FROM freq
     )
     SELECT Gene, Position, Clade, AminoAcid, Freq
     FROM ranked
     WHERE rn = 1 AND Freq >= ?",
    list(subtype, var_type, gene, group_by, clade1, clade2, as.numeric(min_freq))
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)

  c1_dom <- res %>%
    filter(Clade == clade1) %>%
    dplyr::select(Gene, Position, Clade1_AA = AminoAcid, Clade1_Freq = Freq)
  c2_dom <- res %>%
    filter(Clade == clade2) %>%
    dplyr::select(Gene, Position, Clade2_AA = AminoAcid, Clade2_Freq = Freq)

  inner_join(c1_dom, c2_dom, by = c("Gene", "Position")) %>%
    filter(Clade1_AA != Clade2_AA)
}

usage_position_distribution <- function(subtype, var_type, gene, group_by, position, hide_empty_years = FALSE) {
  res <- usage_pairwise_gene_data(subtype, var_type, gene, group_by)
  if (is.null(res)) return(NULL)

  res <- res %>%
    filter(Position == position, !(AminoAcid %in% c("X", "-")))

  if (is.null(res) || nrow(res) == 0) return(NULL)

  has_codon <- "Codon_Usage" %in% colnames(res)
  out <- res %>%
    group_by(Clade, AminoAcid) %>%
    summarise(
      Count = sum(Count, na.rm = TRUE),
      Codon_Usage = if (has_codon) dplyr::first(Codon_Usage) else NA_character_,
      .groups = "drop_last"
    ) %>%
    mutate(Total_in_Clade = sum(Count)) %>%
    mutate(`Frequency(%)` = (Count / Total_in_Clade) * 100) %>%
    ungroup()

  if (group_by == "Year" && isTRUE(hide_empty_years)) {
    out <- out %>% filter(Total_in_Clade > 0)
  }
  if (all(is.na(out$Codon_Usage))) out$Codon_Usage <- NULL
  out
}

usage_entropy_data <- function(subtype, var_type, gene, group_by, clade = "All") {
  clade_filter <- ""
  params <- list(subtype, var_type, gene, group_by)
  if (!identical(clade, "All")) {
    clade_filter <- " AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END = ?"
    params <- c(params, list(clade))
  }

  res <- usage_query(
    paste0(
      "SELECT Position, AminoAcid, SUM(Count) AS AA_Sum
       FROM usage
       WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?",
      clade_filter,
      " AND AminoAcid NOT IN ('X', '-')
       GROUP BY Position, AminoAcid"
    ),
    params
  )
  if (is.null(res)) return(NULL)

  res %>%
    group_by(Position) %>%
    mutate(Pos_Total = sum(AA_Sum), p = AA_Sum / Pos_Total) %>%
    filter(p > 0) %>%
    summarise(
      Entropy = -sum(p * log2(p)),
      Pos_Total = first(Pos_Total),
      .groups = "drop"
    )
}

usage_lollipop_consensus <- function(subtype, var_type, gene, group_by, ref_group, tar_group, min_freq) {
  res <- usage_query(
    "WITH agg AS (
       SELECT Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END AS Clade,
         AminoAcid, SUM(Count) AS Count
       FROM usage
       WHERE \"Group\" = ? AND Variation_Type = ? AND Gene = ? AND Grouping_Type = ?
         AND CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END IN (?, ?)
         AND AminoAcid NOT IN ('X', '-')
       GROUP BY Position,
         CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month_Filter ELSE Clade END,
         AminoAcid
     ),
     freq AS (
       SELECT *,
         SUM(Count) OVER (PARTITION BY Position, Clade) AS Valid_Total,
         100.0 * Count / SUM(Count) OVER (PARTITION BY Position, Clade) AS New_Frequency
       FROM agg
     ),
     ranked AS (
       SELECT *,
         ROW_NUMBER() OVER (PARTITION BY Position, Clade ORDER BY New_Frequency DESC, AminoAcid) AS rn
       FROM freq
     )
     SELECT Position, Clade, AminoAcid, New_Frequency
     FROM ranked
     WHERE rn = 1 AND New_Frequency >= ?",
    list(subtype, var_type, gene, group_by, ref_group, tar_group, as.numeric(min_freq))
  )
  if (is.null(res)) return(NULL)
  res
}

# ==========================================
# 5. MULTI-PATHOGEN ADAPTERS
# ==========================================
# Keep FLU behavior delegated to the original helper functions. Non-FLU
# subtypes are encoded as PATHOGEN:SUBTYPE so existing server call sites can
# route through the same helper names without changing FLU values.
flu_usage_available_genes <- usage_available_genes
flu_usage_available_groups <- usage_available_groups
flu_usage_distinct_group_values <- usage_distinct_group_values
flu_usage_max_position <- usage_max_position
flu_usage_year_month_choices <- usage_year_month_choices
flu_usage_single_position <- usage_single_position
flu_usage_pairwise_gene_data <- usage_pairwise_gene_data
flu_usage_pairwise_differences_for_gene <- usage_pairwise_differences_for_gene
flu_usage_position_distribution <- usage_position_distribution
flu_usage_entropy_data <- usage_entropy_data
flu_usage_lollipop_consensus <- usage_lollipop_consensus
flu_usage_duckdb_available <- usage_duckdb_available
flu_usage_query <- usage_query
flu_usage_sql_in_values <- usage_sql_in_values

PATHOGEN_ADAPTERS <- list(
  FLU = list(
    id = "FLU",
    label = "Influenza",
    schema = "flu",
    subtype_choices = stats::setNames(sort_flu_subtypes(unique(c(metadata_groups, SUBTYPES))), sort_flu_subtypes(unique(c(metadata_groups, SUBTYPES))))
  ),
  RSV = list(
    id = "RSV",
    label = "RSV",
    schema = "standard",
    duckdb = file.path("data", "cache", "RSV", "rsv_explorer.duckdb"),
    metadata = file.path("data", "cache", "RSV", "metadata_global.rds"),
    subtype_choices = c("A" = "RSV:A", "B" = "RSV:B"),
    default_grouping = "Clade"
  ),
  COVID = list(
    id = "COVID",
    label = "COVID",
    schema = "covid",
    duckdb = file.path("data", "cache", "COVID", "covid_explorer.duckdb"),
    metadata = file.path("data", "cache", "COVID", "metadata_summary.rds"),
    subtype_choices = c("SARS-CoV-2" = "COVID:SARS-CoV-2"),
    default_grouping = "Nextstrain_clade"
  ),
  CHIKV = list(
    id = "CHIKV",
    label = "CHIKV",
    schema = "standard",
    duckdb = file.path("data", "cache", "CHIKV", "aa_explorer.duckdb"),
    metadata = file.path("data", "cache", "CHIKV", "metadata_global.rds"),
    subtype_choices = c("CHIKV" = "CHIKV:CHIKV"),
    default_grouping = "Clade"
  )
)

ADAPTER_DISPLAY_GROUP_LIMITS <- c(Nextclade_pango = 100L)

pathogen_choices <- function() {
  stats::setNames(names(PATHOGEN_ADAPTERS), vapply(PATHOGEN_ADAPTERS, `[[`, character(1), "label"))
}

pathogen_subtype_choices <- function(pathogen_id = "FLU") {
  if (is.null(pathogen_id) || length(pathogen_id) == 0 || is.na(pathogen_id[[1]]) || identical(pathogen_id[[1]], "")) {
    pathogen_id <- "FLU"
  }
  cfg <- PATHOGEN_ADAPTERS[[pathogen_id]]
  if (is.null(cfg)) cfg <- PATHOGEN_ADAPTERS$FLU
  cfg$subtype_choices
}

dataset_insights_cache_path <- function(pathogen_id) {
  file.path("data", "cache", pathogen_id, "dataset_insights.rds")
}

format_sequence_count <- function(value) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) return("0")
  if (is.character(value)) return(value[[1]])
  format(as.numeric(value[[1]]), big.mark = ",", scientific = FALSE)
}

clean_year_values <- function(values) {
  suppressWarnings(as.numeric(as.character(values)))
}

unknown_like <- function(values) {
  is_unknown_metadata_value(values)
}

insights_from_summary <- function(pathogen_id, total_sequences, countries_represented, time_range, metadata_groups, metadata_grouping_cols, summary_stats) {
  if (!"YearMonth" %in% names(summary_stats)) summary_stats$YearMonth <- NA_character_
  if (!"region" %in% names(summary_stats)) summary_stats$region <- "Unknown"
  if (!"country" %in% names(summary_stats)) summary_stats$country <- "Unknown"
  if (!"n" %in% names(summary_stats)) summary_stats$n <- 1L

  summary_stats$Year <- clean_year_values(summary_stats$Year)
  summary_stats$Group <- as.character(summary_stats$Group)

  time_plot <- summary_stats %>%
    filter(!is.na(.data$Year)) %>%
    group_by(.data$Year, .data$Group) %>%
    summarise(Count = sum(.data$n, na.rm = TRUE), .groups = "drop")

  geo_plot <- summary_stats %>%
    filter(!is.na(.data$Year), !unknown_like(.data$region)) %>%
    group_by(.data$Year, region = .data$region) %>%
    summarise(Count = sum(.data$n, na.rm = TRUE), .groups = "drop")

  fill_cols <- unique(c(metadata_grouping_cols, "region", "country"))
  fill_cols <- fill_cols[fill_cols %in% names(summary_stats)]
  breakdowns <- list()
  for (time_col in c("Year", "YearMonth")) {
    if (!time_col %in% names(summary_stats)) next
    time_df <- summary_stats
    if (identical(time_col, "Year")) {
      time_df <- time_df %>% filter(!is.na(.data$Year))
    } else {
      time_df <- time_df %>% filter(!is.na(.data$YearMonth), .data$YearMonth != "")
    }
    for (fill_col in fill_cols) {
      df <- time_df %>%
        transmute(
          Group = .data$Group,
          XValue = as.character(.data[[time_col]]),
          FillValue = ifelse(unknown_like(.data[[fill_col]]), "Unknown", as.character(.data[[fill_col]])),
          Count = .data$n
        ) %>%
        group_by(.data$Group, .data$XValue, .data$FillValue) %>%
        summarise(Count = sum(.data$Count, na.rm = TRUE), .groups = "drop") %>%
        mutate(XOrder = .data$XValue)
      breakdowns[[paste(time_col, fill_col, sep = "__")]] <- df
    }
  }

  list(
    pathogen_id = pathogen_id,
    total_sequences = total_sequences,
    countries_represented = countries_represented,
    time_range = time_range,
    metadata_groups = metadata_groups,
    metadata_grouping_cols = metadata_grouping_cols,
    metadata_summary_stats = summary_stats,
    time_plot = time_plot,
    geo_plot = geo_plot,
    breakdowns = breakdowns
  )
}

build_flu_dataset_insights <- function() {
  insights_from_summary(
    pathogen_id = "FLU",
    total_sequences = total_raw,
    countries_represented = total_countries_val,
    time_range = time_range_val,
    metadata_groups = metadata_groups,
    metadata_grouping_cols = metadata_grouping_cols,
    summary_stats = metadata_summary_stats
  )
}

metadata_global_to_summary <- function(pathogen_id, metadata_cache_path, subtype_choices) {
  cache <- readRDS(metadata_cache_path)
  metadata <- cache$metadata_global
  subtype_map <- stats::setNames(unname(subtype_choices), names(subtype_choices))

  if (!"Group" %in% names(metadata)) metadata$Group <- names(subtype_map)[[1]]
  metadata$Group <- as.character(metadata$Group)
  metadata$Group <- ifelse(metadata$Group %in% names(subtype_map), subtype_map[metadata$Group], paste(pathogen_id, metadata$Group, sep = ":"))

  if (!"Year" %in% names(metadata)) metadata$Year <- substr(as.character(metadata$date %||% NA_character_), 1, 4)
  metadata$Year <- clean_year_values(metadata$Year)
  if (!"YearMonth" %in% names(metadata)) {
    if ("YM" %in% names(metadata)) metadata$YearMonth <- metadata$YM
    else if ("date" %in% names(metadata)) metadata$YearMonth <- substr(as.character(metadata$date), 1, 7)
    else metadata$YearMonth <- NA_character_
  }
  if (!"region" %in% names(metadata)) metadata$region <- "Unknown"
  if (!"country" %in% names(metadata)) metadata$country <- "Unknown"

  grouping_cols <- intersect(c("clade", "G_clade", "group_1", "group_2", "group_3", "group_4"), names(metadata))
  for (col in grouping_cols) metadata[[col]] <- as.character(metadata[[col]])

  keep_cols <- unique(c("Year", "YearMonth", "Group", "region", "country", grouping_cols))
  metadata %>%
    dplyr::select(all_of(keep_cols)) %>%
    group_by(across(all_of(keep_cols))) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    list(summary = ., grouping_cols = grouping_cols, cache = cache)
}

build_standard_dataset_insights <- function(pathogen_id) {
  cfg <- PATHOGEN_ADAPTERS[[pathogen_id]]
  info <- metadata_global_to_summary(pathogen_id, cfg$metadata, cfg$subtype_choices)
  summary_stats <- info$summary
  years <- clean_year_values(summary_stats$Year)
  countries <- if ("country" %in% names(summary_stats)) {
    length(unique(summary_stats$country[!unknown_like(summary_stats$country)]))
  } else {
    NA_integer_
  }

  insights_from_summary(
    pathogen_id = pathogen_id,
    total_sequences = info$cache$total_parsed %||% format_sequence_count(sum(summary_stats$n, na.rm = TRUE)),
    countries_represented = countries,
    time_range = if (any(!is.na(years))) paste(min(years, na.rm = TRUE), "-", max(years, na.rm = TRUE)) else "Unavailable",
    metadata_groups = unname(cfg$subtype_choices),
    metadata_grouping_cols = info$grouping_cols,
    summary_stats = summary_stats
  )
}

build_covid_dataset_insights <- function() {
  cfg <- PATHOGEN_ADAPTERS$COVID
  cache <- readRDS(cfg$metadata)
  group_value <- unname(cfg$subtype_choices)[1]
  breakdowns <- lapply(cache$breakdowns, function(df) {
    df <- as.data.frame(df)
    df$Group <- group_value
    df
  })

  time_plot <- cache$time_plot %>%
    transmute(Year = clean_year_values(.data$Year), Group = as.character(.data$clade_who), Count = as.numeric(.data$Count))

  geo_plot <- if ("Year__region" %in% names(breakdowns)) {
    breakdowns$Year__region %>%
      transmute(Year = clean_year_values(.data$XValue), region = as.character(.data$FillValue), Count = as.numeric(.data$Count))
  } else {
    cache$region_plot %>%
      transmute(Year = NA_real_, region = as.character(.data$region), Count = as.numeric(.data$Count))
  }

  grouping_cols <- c("clade_who", "Nextclade_pango", "Nextstrain_clade", "pango_lineage", "region", "country", "division", "host")
  grouping_cols <- grouping_cols[vapply(paste("Year", grouping_cols, sep = "__"), function(key) key %in% names(breakdowns), logical(1))]

  list(
    pathogen_id = "COVID",
    total_sequences = format_sequence_count(cache$global_summary$total_sequences),
    countries_represented = cache$global_summary$countries_represented,
    time_range = cache$global_summary$time_span,
    metadata_groups = group_value,
    metadata_grouping_cols = grouping_cols,
    metadata_summary_stats = data.frame(),
    time_plot = time_plot,
    geo_plot = geo_plot,
    breakdowns = breakdowns
  )
}

build_dataset_insights_cache <- function(pathogen_id) {
  insights <- switch(
    pathogen_id,
    FLU = build_flu_dataset_insights(),
    COVID = build_covid_dataset_insights(),
    RSV = build_standard_dataset_insights("RSV"),
    CHIKV = build_standard_dataset_insights("CHIKV"),
    stop("Unsupported pathogen for Dataset Insights: ", pathogen_id, call. = FALSE)
  )
  path <- dataset_insights_cache_path(pathogen_id)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(insights, path)
  insights
}

ensure_dataset_insights_cache <- function(pathogen_id) {
  path <- dataset_insights_cache_path(pathogen_id)
  if (file.exists(path)) return(path)
  build_dataset_insights_cache(pathogen_id)
  path
}

ensure_all_dataset_insights_cache <- function() {
  invisible(vapply(unname(pathogen_choices()), ensure_dataset_insights_cache, character(1)))
}

load_dataset_insights <- function(pathogen_id) {
  path <- ensure_dataset_insights_cache(pathogen_id)
  readRDS(path)
}

pathogen_from_subtype <- function(subtype) {
  subtype <- as.character(subtype)
  if (length(subtype) == 0 || is.na(subtype[[1]]) || identical(subtype[[1]], "")) return(NA_character_)
  if (!grepl(":", subtype[[1]], fixed = TRUE)) return("FLU")
  strsplit(subtype[[1]], ":", fixed = TRUE)[[1]][[1]]
}

adapter_subtype_value <- function(subtype) {
  subtype <- as.character(subtype)
  if (length(subtype) == 0 || is.na(subtype[[1]])) return(subtype)
  parts <- strsplit(subtype[[1]], ":", fixed = TRUE)[[1]]
  if (length(parts) < 2) subtype[[1]] else paste(parts[-1], collapse = ":")
}

adapter_config <- function(subtype) {
  pathogen_id <- pathogen_from_subtype(subtype)
  if (length(pathogen_id) == 0 || is.na(pathogen_id) || !pathogen_id %in% names(PATHOGEN_ADAPTERS)) return(NULL)
  PATHOGEN_ADAPTERS[[pathogen_id]]
}

is_flu_subtype <- function(subtype) {
  identical(pathogen_from_subtype(subtype), "FLU")
}

missing_subtype <- function(subtype) {
  subtype <- as.character(subtype)
  length(subtype) == 0 || is.na(subtype[[1]]) || identical(subtype[[1]], "")
}

adapter_db_env <- new.env(parent = emptyenv())

adapter_db_conn <- function(cfg) {
  if (!isTRUE(USE_DUCKDB) || is.null(cfg$duckdb) || !file.exists(cfg$duckdb)) return(NULL)
  key <- cfg$id
  con <- adapter_db_env[[key]]
  if (is.null(con) || !DBI::dbIsValid(con)) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = cfg$duckdb, read_only = TRUE)
    DBI::dbExecute(con, "PRAGMA memory_limit='700MB'")
    adapter_db_env[[key]] <- con
  }
  con
}

adapter_query <- function(cfg, sql, params = NULL) {
  con <- adapter_db_conn(cfg)
  if (is.null(con)) return(NULL)
  DBI::dbGetQuery(con, sql, params = params)
}

adapter_quote <- function(cfg, value) {
  con <- adapter_db_conn(cfg)
  if (is.null(con)) return("NULL")
  as.character(DBI::dbQuoteString(con, value))
}

adapter_sql_in_values <- function(cfg, values) {
  con <- adapter_db_conn(cfg)
  if (is.null(con) || length(values) == 0) return(NULL)
  paste(DBI::dbQuoteString(con, values), collapse = ", ")
}

standard_group_value_sql <- function() {
  "CASE WHEN Grouping_Type = 'Year_Month' THEN Year_Month ELSE Clade END"
}

adapter_available_genes <- function(subtype, var_type) {
  if (!identical(var_type, "AA")) return(character(0))
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(character(0))
  raw_subtype <- adapter_subtype_value(subtype)
  if (identical(cfg$schema, "covid")) {
    res <- adapter_query(cfg, "SELECT DISTINCT Protein AS Gene FROM usage ORDER BY Gene")
  } else {
    res <- adapter_query(cfg, "SELECT DISTINCT Gene FROM usage WHERE \"Group\" = ? ORDER BY Gene", list(raw_subtype))
  }
  if (is.null(res) || nrow(res) == 0) character(0) else sort(stats::na.omit(as.character(res$Gene)))
}

adapter_available_groups <- function(subtype, var_type, gene) {
  if (!identical(var_type, "AA")) return(character(0))
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(character(0))
  raw_subtype <- adapter_subtype_value(subtype)
  if (identical(cfg$schema, "covid")) {
    res <- adapter_query(cfg, "SELECT DISTINCT GroupType AS Grouping_Type FROM usage WHERE Protein = ? ORDER BY Grouping_Type", list(gene))
  } else {
    res <- adapter_query(cfg, "SELECT DISTINCT Grouping_Type FROM usage WHERE \"Group\" = ? AND Gene = ? ORDER BY Grouping_Type", list(raw_subtype, gene))
  }
  if (is.null(res) || nrow(res) == 0) return(character(0))
  groups <- sort(stats::na.omit(as.character(res$Grouping_Type)))
  preferred <- cfg$default_grouping
  if (!is.null(preferred) && preferred %in% groups) c(preferred, setdiff(groups, preferred)) else groups
}

adapter_distinct_group_values <- function(subtype, var_type, gene, group_by) {
  if (!identical(var_type, "AA")) return(character(0))
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(character(0))
  raw_subtype <- adapter_subtype_value(subtype)
  if (identical(cfg$schema, "covid")) {
    limit <- if (length(group_by) == 1 && group_by %in% names(ADAPTER_DISPLAY_GROUP_LIMITS)) {
      ADAPTER_DISPLAY_GROUP_LIMITS[[group_by]]
    } else {
      NULL
    }
    if (!is.null(limit) && !is.na(limit)) {
      res <- adapter_query(
        cfg,
        "SELECT GroupValue AS Clade, SUM(Count) AS TotalCount
         FROM usage
         WHERE Protein = ? AND GroupType = ? AND GroupValue IS NOT NULL AND GroupValue <> ''
         GROUP BY GroupValue
         ORDER BY TotalCount DESC, GroupValue
         LIMIT ?",
        list(gene, group_by, as.integer(limit))
      )
      if (is.null(res) || nrow(res) == 0) return(character(0))
      return(sort(stats::na.omit(as.character(res$Clade))))
    }
    res <- adapter_query(
      cfg,
      "SELECT DISTINCT GroupValue AS Clade FROM usage WHERE Protein = ? AND GroupType = ? ORDER BY Clade",
      list(gene, group_by)
    )
  } else {
    res <- adapter_query(
      cfg,
      paste0("SELECT DISTINCT ", standard_group_value_sql(), " AS Clade FROM usage WHERE \"Group\" = ? AND Gene = ? AND Grouping_Type = ? ORDER BY Clade"),
      list(raw_subtype, gene, group_by)
    )
  }
  if (is.null(res) || nrow(res) == 0) return(character(0))
  values <- sort(stats::na.omit(as.character(res$Clade)))
  special_values <- c("Unknown", "unassigned", "Unassigned")
  c(setdiff(values, special_values), intersect(special_values, values))
}

adapter_position_choices <- function(subtype, var_type, gene) {
  if (!identical(var_type, "AA")) return(character(0))
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(character(0))
  raw_subtype <- adapter_subtype_value(subtype)
  if (identical(cfg$schema, "covid")) {
    res <- adapter_query(
      cfg,
      "SELECT DISTINCT Position AS Position_Key, COALESCE(PositionLabel, CAST(PositionBase AS VARCHAR)) AS Position_Label,
              PositionBase, InsertionOffset
       FROM usage
       WHERE Protein = ?
       ORDER BY PositionBase, InsertionOffset, Position_Label",
      list(gene)
    )
  } else {
    res <- adapter_query(
      cfg,
      "SELECT DISTINCT Position_Key, Position_Label, Position, Position_Offset
       FROM usage
       WHERE \"Group\" = ? AND Gene = ?
       ORDER BY Position, Position_Offset, Position_Label",
      list(raw_subtype, gene)
    )
  }
  if (is.null(res) || nrow(res) == 0) return(character(0))
  stats::setNames(as.character(res$Position_Key), as.character(res$Position_Label))
}

adapter_position_filter <- function(cfg, position) {
  position_chr <- as.character(position)
  if (identical(cfg$schema, "covid")) {
    if (!is.na(suppressWarnings(as.numeric(position_chr)))) {
      return(list(sql = "PositionBase = ?", value = suppressWarnings(as.numeric(position_chr))))
    }
    return(list(sql = "Position = ?", value = position_chr))
  }
  if (grepl("^\\d{9}__\\d{3}$", position_chr)) {
    return(list(sql = "Position_Key = ?", value = position_chr))
  }
  list(sql = "Position = ?", value = suppressWarnings(as.numeric(position_chr)))
}

adapter_max_position <- function(subtype, var_type, gene) {
  if (!identical(var_type, "AA")) return(NA_real_)
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(NA_real_)
  raw_subtype <- adapter_subtype_value(subtype)
  if (identical(cfg$schema, "covid")) {
    res <- adapter_query(cfg, "SELECT MAX(PositionBase) AS max_position FROM usage WHERE Protein = ?", list(gene))
  } else {
    res <- adapter_query(cfg, "SELECT MAX(Position) AS max_position FROM usage WHERE \"Group\" = ? AND Gene = ?", list(raw_subtype, gene))
  }
  if (is.null(res) || nrow(res) == 0 || is.na(res$max_position[[1]])) NA_real_ else as.numeric(res$max_position[[1]])
}

adapter_year_month_choices <- function(subtype, var_type, gene, group_by, position) {
  if (!identical(var_type, "AA")) return(character(0))
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(character(0))
  raw_subtype <- adapter_subtype_value(subtype)
  if (identical(cfg$schema, "covid")) {
    res <- adapter_query(
      cfg,
      "SELECT DISTINCT GroupValue AS Year_Month FROM usage WHERE Protein = ? AND GroupType = ? AND PositionBase = ? ORDER BY Year_Month",
      list(gene, group_by, as.numeric(position))
    )
  } else {
    res <- adapter_query(
      cfg,
      "SELECT DISTINCT Year_Month FROM usage WHERE \"Group\" = ? AND Gene = ? AND Grouping_Type = ? AND Position = ? AND Year_Month IS NOT NULL ORDER BY Year_Month",
      list(raw_subtype, gene, group_by, as.numeric(position))
    )
  }
  if (is.null(res) || nrow(res) == 0) return(character(0))
  values <- stats::na.omit(as.character(res$Year_Month))
  special_values <- c("Unknown", "unassigned", "Unassigned")
  c(intersect(special_values, values), sort(setdiff(values, special_values)))
}

adapter_pairwise_gene_data <- function(subtype, var_type, gene, group_by, clades = NULL) {
  if (!identical(var_type, "AA")) return(data.frame())
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(data.frame())
  raw_subtype <- adapter_subtype_value(subtype)
  clade_filter <- ""
  params <- list()
  limited_group <- length(group_by) == 1 && group_by %in% names(ADAPTER_DISPLAY_GROUP_LIMITS)
  if (is.null(clades) && identical(cfg$schema, "covid") && isTRUE(limited_group)) {
    clades <- adapter_distinct_group_values(subtype, var_type, gene, group_by)
  }
  if (!is.null(clades) && length(clades) > 0) {
    in_values <- adapter_sql_in_values(cfg, clades)
    value_expr <- if (identical(cfg$schema, "covid")) "GroupValue" else standard_group_value_sql()
    if (!is.null(in_values)) clade_filter <- paste0(" AND ", value_expr, " IN (", in_values, ")")
  }
  group_literal <- adapter_quote(cfg, subtype)

  if (identical(cfg$schema, "covid")) {
    params <- list(gene, group_by)
    sql <- paste0(
      "SELECT ", group_literal, " AS \"Group\", Protein AS Gene, GroupValue AS Clade,
              PositionBase AS Position, AA AS AminoAcid, SUM(Count) AS Count,
              ANY_VALUE(Codon) AS Codon_Usage
       FROM usage
       WHERE Protein = ? AND GroupType = ?",
      clade_filter,
      " GROUP BY Protein, GroupValue, PositionBase, AA"
    )
  } else {
    params <- list(raw_subtype, gene, group_by)
    sql <- paste0(
      "SELECT ", group_literal, " AS \"Group\", Gene, ", standard_group_value_sql(), " AS Clade,
              Position, AminoAcid, SUM(Count) AS Count,
              ANY_VALUE(Codon_Usage) AS Codon_Usage
       FROM usage
       WHERE \"Group\" = ? AND Gene = ? AND Grouping_Type = ?",
      clade_filter,
      " GROUP BY Gene, ", standard_group_value_sql(), ", Position, AminoAcid"
    )
  }
  res <- adapter_query(cfg, sql, params)
  if (is.null(res)) return(NULL)
  if ("Codon_Usage" %in% names(res) && all(is.na(res$Codon_Usage))) res$Codon_Usage <- NULL
  res
}

adapter_single_position <- function(subtype, var_type, gene, group_by, position, allowed_yms = NULL, min_seqs = 1, hide_empty_years = FALSE) {
  cfg <- adapter_config(subtype)
  if (is.null(cfg)) return(data.frame())
  raw_subtype <- adapter_subtype_value(subtype)
  pos_filter <- adapter_position_filter(cfg, position)
  if (identical(cfg$schema, "covid")) {
    data <- adapter_query(
      cfg,
      paste0(
        "SELECT ", adapter_quote(cfg, subtype), " AS \"Group\", Protein AS Gene, GroupValue AS Clade,
                PositionBase AS Position, COALESCE(PositionLabel, CAST(PositionBase AS VARCHAR)) AS Position_Label,
                Position AS Position_Key, AA AS AminoAcid, SUM(Count) AS Count, ANY_VALUE(Codon) AS Codon_Usage
         FROM usage
         WHERE Protein = ? AND GroupType = ? AND ", pos_filter$sql, " AND AA NOT IN ('X', '-')
         GROUP BY Protein, GroupValue, PositionBase, PositionLabel, Position, AA"
      ),
      list(gene, group_by, pos_filter$value)
    )
  } else {
    data <- adapter_query(
      cfg,
      paste0(
        "SELECT ", adapter_quote(cfg, subtype), " AS \"Group\", Gene, ", standard_group_value_sql(), " AS Clade,
                Position, Position_Label, Position_Key, AminoAcid, SUM(Count) AS Count, ANY_VALUE(Codon_Usage) AS Codon_Usage
         FROM usage
         WHERE \"Group\" = ? AND Gene = ? AND Grouping_Type = ? AND ", pos_filter$sql, " AND AminoAcid NOT IN ('X', '-')
         GROUP BY Gene, ", standard_group_value_sql(), ", Position, Position_Label, Position_Key, AminoAcid"
      ),
      list(raw_subtype, gene, group_by, pos_filter$value)
    )
  }
  if (is.null(data)) return(NULL)
  if (nrow(data) == 0) return(data.frame())
  if (!is.null(allowed_yms) && length(allowed_yms) > 0 && group_by == "Year_Month") {
    data <- data %>% filter(.data$Clade %in% allowed_yms)
  }
  if (nrow(data) == 0) {
    out <- data.frame(Group=character(), Gene=character(), Position=numeric(), AminoAcid=character(), Count=numeric(), Valid_Total=numeric(), `Frequency(%)`=numeric(), check.names = FALSE)
    out[[group_by]] <- character()
    return(out)
  }
  out <- data %>%
    group_by(.data$Group, .data$Gene, .data$Position, Clade = .data$Clade, .data$AminoAcid) %>%
    summarise(Count = sum(.data$Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data$Clade) %>%
    mutate(Valid_Total = sum(.data$Count, na.rm = TRUE), `Frequency(%)` = (.data$Count / .data$Valid_Total) * 100) %>%
    ungroup() %>%
    filter(.data$Valid_Total >= min_seqs) %>%
    rename(!!group_by := Clade)
  if (group_by == "Year" && isTRUE(hide_empty_years)) out <- out %>% filter(.data$Valid_Total > 0)
  out
}

adapter_pairwise_differences_for_gene <- function(subtype, var_type, gene, group_by, clade1, clade2, min_freq) {
  res <- adapter_pairwise_gene_data(subtype, var_type, gene, group_by, clades = c(clade1, clade2))
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res <- res %>% filter(!(.data$AminoAcid %in% c("X", "-")))
  top <- res %>%
    group_by(.data$Gene, .data$Position, .data$Clade, .data$AminoAcid) %>%
    summarise(Variant_Count = sum(.data$Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data$Gene, .data$Position, .data$Clade) %>%
    mutate(Total_Seqs = sum(.data$Variant_Count, na.rm = TRUE), Freq = 100 * .data$Variant_Count / .data$Total_Seqs) %>%
    arrange(.data$Gene, .data$Position, .data$Clade, desc(.data$Freq), .data$AminoAcid) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    filter(.data$Freq >= min_freq)
  if (nrow(top) == 0) return(NULL)
  c1_dom <- top %>% filter(.data$Clade == clade1) %>% dplyr::select(Gene, Position, Clade1_AA = AminoAcid, Clade1_Freq = Freq)
  c2_dom <- top %>% filter(.data$Clade == clade2) %>% dplyr::select(Gene, Position, Clade2_AA = AminoAcid, Clade2_Freq = Freq)
  inner_join(c1_dom, c2_dom, by = c("Gene", "Position")) %>% filter(.data$Clade1_AA != .data$Clade2_AA)
}

adapter_position_distribution <- function(subtype, var_type, gene, group_by, position, hide_empty_years = FALSE) {
  res <- adapter_single_position(subtype, var_type, gene, group_by, position, min_seqs = 1, hide_empty_years = hide_empty_years)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  group_col <- group_by
  if (!group_col %in% names(res)) return(NULL)
  res <- res %>% dplyr::rename(Clade = !!sym(group_col))
  if (nrow(res) == 0) return(NULL)
  out <- res %>%
    group_by(.data$Clade, .data$AminoAcid) %>%
    summarise(Count = sum(.data$Count, na.rm = TRUE), .groups = "drop_last") %>%
    mutate(Total_in_Clade = sum(.data$Count), `Frequency(%)` = (.data$Count / .data$Total_in_Clade) * 100) %>%
    ungroup()
  if (group_by == "Year" && isTRUE(hide_empty_years)) out <- out %>% filter(.data$Total_in_Clade > 0)
  out
}

adapter_entropy_data <- function(subtype, var_type, gene, group_by, clade = "All") {
  res <- adapter_pairwise_gene_data(subtype, var_type, gene, group_by, clades = if (identical(clade, "All")) NULL else clade)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res %>%
    filter(!(.data$AminoAcid %in% c("X", "-"))) %>%
    group_by(.data$Position, .data$AminoAcid) %>%
    summarise(AA_Sum = sum(.data$Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data$Position) %>%
    mutate(Pos_Total = sum(.data$AA_Sum), p = .data$AA_Sum / .data$Pos_Total) %>%
    filter(.data$p > 0) %>%
    summarise(
      Entropy = -sum(.data$p * log2(.data$p)),
      Pos_Total = first(.data$Pos_Total),
      .groups = "drop"
    )
}

usage_available_genes <- function(subtype, var_type) {
  if (missing_subtype(subtype)) return(character(0))
  if (is_flu_subtype(subtype)) flu_usage_available_genes(subtype, var_type) else adapter_available_genes(subtype, var_type)
}

usage_available_groups <- function(subtype, var_type, gene) {
  if (missing_subtype(subtype)) return(character(0))
  if (is_flu_subtype(subtype)) flu_usage_available_groups(subtype, var_type, gene) else adapter_available_groups(subtype, var_type, gene)
}

usage_distinct_group_values <- function(subtype, var_type, gene, group_by) {
  if (missing_subtype(subtype)) return(character(0))
  if (is_flu_subtype(subtype)) flu_usage_distinct_group_values(subtype, var_type, gene, group_by) else adapter_distinct_group_values(subtype, var_type, gene, group_by)
}

usage_max_position <- function(subtype, var_type, gene) {
  if (missing_subtype(subtype)) return(NA_real_)
  if (is_flu_subtype(subtype)) flu_usage_max_position(subtype, var_type, gene) else adapter_max_position(subtype, var_type, gene)
}

usage_year_month_choices <- function(subtype, var_type, gene, group_by, position) {
  if (missing_subtype(subtype)) return(character(0))
  if (is_flu_subtype(subtype)) flu_usage_year_month_choices(subtype, var_type, gene, group_by, position) else adapter_year_month_choices(subtype, var_type, gene, group_by, position)
}

usage_position_choices <- function(subtype, var_type, gene) {
  if (missing_subtype(subtype)) return(character(0))
  if (is_flu_subtype(subtype)) flu_position_choices(subtype, var_type, gene) else adapter_position_choices(subtype, var_type, gene)
}

usage_single_position <- function(subtype, var_type, gene, group_by, position, allowed_yms = NULL, min_seqs = 1, hide_empty_years = FALSE) {
  if (missing_subtype(subtype)) return(data.frame())
  if (is_flu_subtype(subtype)) {
    flu_usage_single_position(subtype, var_type, gene, group_by, position, allowed_yms, min_seqs, hide_empty_years)
  } else {
    adapter_single_position(subtype, var_type, gene, group_by, position, allowed_yms, min_seqs, hide_empty_years)
  }
}

usage_pairwise_gene_data <- function(subtype, var_type, gene, group_by, clades = NULL) {
  if (missing_subtype(subtype)) return(data.frame())
  if (is_flu_subtype(subtype)) flu_usage_pairwise_gene_data(subtype, var_type, gene, group_by, clades) else adapter_pairwise_gene_data(subtype, var_type, gene, group_by, clades)
}

usage_pairwise_differences_for_gene <- function(subtype, var_type, gene, group_by, clade1, clade2, min_freq) {
  if (missing_subtype(subtype)) return(NULL)
  if (is_flu_subtype(subtype)) {
    flu_usage_pairwise_differences_for_gene(subtype, var_type, gene, group_by, clade1, clade2, min_freq)
  } else {
    adapter_pairwise_differences_for_gene(subtype, var_type, gene, group_by, clade1, clade2, min_freq)
  }
}

usage_position_distribution <- function(subtype, var_type, gene, group_by, position, hide_empty_years = FALSE) {
  if (missing_subtype(subtype)) return(data.frame())
  if (is_flu_subtype(subtype)) flu_usage_position_distribution(subtype, var_type, gene, group_by, position, hide_empty_years) else adapter_position_distribution(subtype, var_type, gene, group_by, position, hide_empty_years)
}

usage_entropy_data <- function(subtype, var_type, gene, group_by, clade = "All") {
  if (missing_subtype(subtype)) return(NULL)
  if (is_flu_subtype(subtype)) flu_usage_entropy_data(subtype, var_type, gene, group_by, clade) else adapter_entropy_data(subtype, var_type, gene, group_by, clade)
}

usage_lollipop_consensus <- function(subtype, var_type, gene, group_by, ref_group, tar_group, min_freq) {
  if (missing_subtype(subtype)) return(NULL)
  if (is_flu_subtype(subtype)) {
    flu_usage_lollipop_consensus(subtype, var_type, gene, group_by, ref_group, tar_group, min_freq)
  } else {
    adapter_pairwise_differences_for_gene(subtype, var_type, gene, group_by, ref_group, tar_group, min_freq)
  }
}
