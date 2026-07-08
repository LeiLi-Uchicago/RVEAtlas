#!/usr/bin/env Rscript
# =============================================================================
# fix_covid_year_month.R
#
# Corrects a double-counting bug in the COVID usage cache
# (data/cache/COVID/covid_explorer.duckdb).
#
# The GroupType='Year_month' rows were aggregated over BOTH full-coverage clade
# schemes (Nextstrain_clade + Nextclade_pango), so every sequence is counted
# TWICE. Result: for any position the Year-Month "total counted AAs" is exactly
# 2x the true value (~18M instead of ~9M), while clade groupings are correct.
# The AA proportions are unaffected (the doubling is perfectly uniform: every
# cell aggregates to an even count), so only absolute counts/totals are wrong.
#
# Fix: replace the Year_month rows with their aggregate at the natural key, with
# Count halved. Verified beforehand that every codon-level cell sum is even, so
# halving is exact and lossless.
#
# Guarded + idempotent: only applies when Year_month total ~= 2x a clade total;
# re-running on already-fixed data is a no-op. NOTE: the covid_explorer.duckdb
# ships prebuilt as release data, so the ROOT cause is in the upstream build
# pipeline (outside this repo) -- fix it there too, or a re-download reintroduces
# this and you must re-run this script.
#
#   Rscript tools/fix_covid_year_month.R
# =============================================================================
suppressWarnings(suppressMessages(library(DBI)))

DB <- "data/cache/COVID/covid_explorer.duckdb"
KEY <- c("Protein", "GroupType", "GroupValue", "Position", "PositionLabel",
         "PositionBase", "InsertionOffset", "PositionOrder", "AA", "Codon")

con <- dbConnect(duckdb::duckdb(), DB, read_only = FALSE)
on.exit(dbDisconnect(con, shutdown = TRUE))

total <- function(gt) dbGetQuery(con, sprintf(
  "SELECT COALESCE(SUM(Count),0) AS x FROM usage WHERE GroupType='%s'", gt))$x

ym <- total("Year_month")
ns <- total("Nextstrain_clade")
cat(sprintf("Before: Year_month total = %.0f ; Nextstrain_clade total = %.0f ; ratio = %.4f\n",
            ym, ns, ym / ns))

if (abs(ym / ns - 1) < 0.01) {
  cat("Year_month already ~= clade total -- looks already fixed. No changes made.\n")
} else if (abs(ym / ns - 2) < 0.01) {
  # Guard: confirm halving is integer-exact at the full (codon-level) key.
  odd <- dbGetQuery(con, sprintf(
    "SELECT COALESCE(SUM(CASE WHEN c %% 2 <> 0 THEN 1 ELSE 0 END),0) AS odd FROM (
       SELECT %s, SUM(Count) c FROM usage WHERE GroupType='Year_month' GROUP BY %s)",
    paste(KEY, collapse = ", "), paste(KEY, collapse = ", ")))$odd
  if (odd > 0) stop(sprintf("Aborting: %d Year_month cells have odd totals; halving would not be exact.", odd))

  cat("Applying fix (aggregate Year_month to natural key, halve Count)...\n")
  dbExecute(con, "BEGIN TRANSACTION")
  ok <- tryCatch({
    dbExecute(con, sprintf(
      "CREATE TEMP TABLE ym_fixed AS
         SELECT %s, SUM(Count) / 2.0 AS Count FROM usage
         WHERE GroupType='Year_month' GROUP BY %s",
      paste(KEY, collapse = ", "), paste(KEY, collapse = ", ")))
    dbExecute(con, "DELETE FROM usage WHERE GroupType='Year_month'")
    dbExecute(con, "INSERT INTO usage SELECT * FROM ym_fixed")
    dbExecute(con, "DROP TABLE ym_fixed")
    dbExecute(con, "COMMIT")
    TRUE
  }, error = function(e) { dbExecute(con, "ROLLBACK"); message("Rolled back: ", conditionMessage(e)); FALSE })
  if (!ok) stop("Migration failed; database unchanged.")

  ym2 <- total("Year_month")
  cat(sprintf("After:  Year_month total = %.0f ; Nextstrain_clade total = %.0f ; ratio = %.4f\n",
              ym2, ns, ym2 / ns))
  if (abs(ym2 / ns - 1) > 0.01) stop("Verification FAILED: Year_month total does not match clade total.")
  cat("Fix applied and verified: Year_month now matches the clade totals.\n")
} else {
  stop(sprintf("Unexpected Year_month/clade ratio (%.4f); not the known 2x doubling. No changes made.", ym / ns))
}
