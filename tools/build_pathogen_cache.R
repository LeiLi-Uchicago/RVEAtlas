#!/usr/bin/env Rscript
# =============================================================================
# build_pathogen_cache.R
#
# Rebuild EVERY cache for one pathogen from raw, in one offline step:
#   duckdb + metadata  ->  conservation_entropy.rds  ->  dataset_insights.rds
#
# Why: rebuilding only the duckdb leaves the pre-computed conservation-entropy
# cache stale (it stores a signature of the source duckdb files). The FIRST app
# launch afterwards then silently recomputes entropy for the pathogen(s) -- ~15s
# for all four -- AFTER "RDS cache loaded successfully." and BEFORE the browser
# opens. Running this build step offline keeps launches fast. This mirrors what
# tools/build_covid_cache.R already does for COVID, now available for every
# pathogen so their rebuild behaviour is consistent.
#
#   Rscript tools/build_pathogen_cache.R <FLU|RSV|COVID|CHIKV>
# =============================================================================

args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript tools/build_pathogen_cache.R <FLU|RSV|COVID|CHIKV>")
}
pid   <- toupper(args[[1]])
valid <- c("FLU", "RSV", "COVID", "CHIKV")
if (!pid %in% valid) {
  stop("Unknown pathogen '", pid, "'. Expected one of: ", paste(valid, collapse = ", "))
}

# COVID already ships a complete standalone builder (duckdb + metadata + aux),
# so delegate to it rather than duplicating the covid-schema logic here.
if (identical(pid, "COVID")) {
  status <- system2("Rscript", c(shQuote(file.path("tools", "build_covid_cache.R")), "all"))
  quit(status = if (identical(as.integer(status), 0L)) 0L else 1L)
}

# FLU rebuilds its duckdb + RDS from raw only when this env var is set, and it
# must be set BEFORE global.R is sourced (the FLU build runs while sourcing).
if (identical(pid, "FLU")) Sys.setenv(FLUEXPLORER_REBUILD_FLU_CACHE = "true")

message("Loading app functions (global.R)...")
suppressMessages(source("global.R", chdir = FALSE))

# Standard-schema pathogens (RSV, CHIKV): force a fresh duckdb + metadata from
# raw. (FLU's duckdb/RDS were already rebuilt above while sourcing global.R.)
if (pid %in% c("RSV", "CHIKV")) {
  # Release the read-only connection global.R's startup opened, so the rebuild
  # can atomically replace the duckdb file.
  if (exists(pid, envir = adapter_db_env, inherits = FALSE)) {
    con <- get(pid, envir = adapter_db_env)
    if (!is.null(con) && DBI::dbIsValid(con)) {
      try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    }
    rm(list = pid, envir = adapter_db_env)
  }
  message("Rebuilding ", pid, " duckdb + metadata from raw...")
  build_standard_adapter_cache(pid)
}

# Pre-compute the derived caches so the app never rebuilds them at launch.
# FLU's conservation was already refreshed by global.R's startup tail (its duckdb
# changed during sourcing), so only force it for the standard pathogens whose
# duckdb we just replaced afterwards.
force_cons <- !identical(pid, "FLU")
message("Rebuilding conservation entropy cache for ", pid, "...")
ensure_conservation_entropy_cache(pid, force = force_cons)
cons <- readRDS(conservation_cache_path(pid))
message("  conservation rows: ", format(nrow(cons$data), big.mark = ","))

message("Rebuilding dataset-insights cache for ", pid, "...")
di_path <- dataset_insights_cache_path(pid)
if (file.exists(di_path)) file.remove(di_path)
di <- load_dataset_insights(pid)
message("  dataset-insights total sequences: ",
        format(di$total_sequences %||% NA, big.mark = ","))

message("Done: ", pid, " caches rebuilt (duckdb + metadata + conservation + dataset-insights).")
