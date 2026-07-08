#!/usr/bin/env Rscript
# =============================================================================
# generate_glyco_epitopes.R
#
# Writes the "N-linked glycosylation" epitope group into the structure epitope
# TSVs for the genes that have a 3D structure:
#   HA  -> www/structures/epitopes_h1.tsv / _h3.tsv / _h5.tsv
#   RSV F -> www/structures/epitopes_rsv_f.tsv   (subtypes RSV:A, RSV:B)
#   SARS-CoV-2 S -> www/structures/epitopes_cov_s.tsv
#
# Sites = potential N-X-S/T sequons (X != Pro) observed across the sequence
# history. A residue counts as "observed" at a position if, in ANY time bucket
# (year, or year-month for COVID) with at least `min_seqs` sequences, its
# within-bucket frequency is >= `q`. Union across buckets keeps sites that were
# common only in one era (e.g. 1968 H3 founder sites, or variant-specific Spike
# glycans) instead of washing them out by pooling.
#
# Numbering: HA maps app Full_HA_Position -> HA_Region + Numbering_Position via
# ha_numbering_review_table.csv. RSV F and SARS-CoV-2 S use identity numbering
# (app position already equals the epitope/structure numbering), so the flagged
# position is written directly under region "F" / "S".
#
# One row per (subtype, region): site "N-Gly", positions comma-separated, matching
# the other epitope groups. Idempotent: existing glyco rows are removed first.
# Re-run whenever a sequence cache is rebuilt.
#
#   Rscript tools/generate_glyco_epitopes.R
# =============================================================================
suppressWarnings(suppressMessages(library(DBI)))

Q         <- 0.20   # within-bucket frequency threshold
MIN_SEQS  <- 20     # ignore time buckets with fewer sequences at a position
GROUP     <- "N-linked glycosylation"
SITE      <- "N-Gly"
COLOR     <- "#009E73"
SOURCE    <- "Potential N-X-S/T sequon across sequence history"
GROUP_COL <- 5L     # group is column 5 in the epitope TSVs

# subtype -> epitope TSV
FILE_FOR <- c(
  H1N1 = "www/structures/epitopes_h1.tsv",
  H3N2 = "www/structures/epitopes_h3.tsv",
  H5NX = "www/structures/epitopes_h5.tsv",
  "RSV:A" = "www/structures/epitopes_rsv_f.tsv",
  "RSV:B" = "www/structures/epitopes_rsv_f.tsv",
  "COVID:SARS-CoV-2" = "www/structures/epitopes_cov_s.tsv"
)

# ---- shared scan core ------------------------------------------------------

# df columns: pos, aa, bucket, n  ->  list(app position -> observed AA vector)
observed_sets <- function(df) {
  df <- df[!is.na(df$aa) & !df$aa %in% c("X", "-", "*") & !is.na(df$bucket), ]
  bt <- tapply(df$n, list(df$pos, df$bucket), sum)
  df$bt <- bt[cbind(as.character(df$pos), as.character(df$bucket))]
  df$freq <- df$n / df$bt
  df <- df[df$bt >= MIN_SEQS & df$freq >= Q, ]
  split(df$aa, df$pos)
}

# Asn app positions starting an N-X-S/T sequon (X != Pro).
scan_positions <- function(sets) {
  pos <- sort(as.integer(names(sets)))
  has <- function(p, aas) { s <- sets[[as.character(p)]]; !is.null(s) && any(aas %in% s) }
  not_pro_only <- function(p) { s <- sets[[as.character(p)]]; !is.null(s) && !all(s == "P") }
  pos[vapply(pos, function(p)
    has(p, "N") && not_pro_only(p + 1) && has(p + 2, c("S", "T")), logical(1))]
}

# ---- per-dataset builders: data.frame(subtype, region, position) -----------

build_flu_ha <- function() {
  con <- dbConnect(duckdb::duckdb(), "data/cache/FLU/flu_explorer.duckdb", read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  num_all <- read.csv("ha_numbering_review_table.csv", stringsAsFactors = FALSE, check.names = FALSE)
  do.call(rbind, lapply(c("H1N1", "H3N2", "H5NX"), function(subtype) {
    d <- dbGetQuery(con, sprintf(
      "SELECT Position pos, AminoAcid aa, Year bucket, SUM(Count) n FROM usage
       WHERE \"Group\"='%s' AND Gene='HA' AND Variation_Type='AA'
       GROUP BY Position, AminoAcid, Year", subtype))
    flagged <- scan_positions(observed_sets(d))
    nm <- num_all[num_all$Subtype == subtype & !as.logical(num_all$Is_Alignment_Gap), ]
    mp <- nm[match(flagged, nm$Full_HA_Position), c("HA_Region", "Numbering_Position")]
    ok <- !is.na(mp$HA_Region) & mp$HA_Region %in% c("HA1", "HA2") & !is.na(mp$Numbering_Position)
    mp <- unique(mp[ok, ])
    data.frame(subtype = subtype, gene = "HA", region = mp$HA_Region,
               position = as.integer(mp$Numbering_Position), stringsAsFactors = FALSE)
  }))
}

build_rsv_f <- function() {
  con <- dbConnect(duckdb::duckdb(), "data/cache/RSV/rsv_explorer.duckdb", read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  # Grouping_Type='Year' avoids multi-counting across the clade/time aggregations.
  do.call(rbind, lapply(c("A", "B"), function(grp) {
    d <- dbGetQuery(con, sprintf(
      "SELECT Position pos, AminoAcid aa, Year bucket, SUM(Count) n FROM usage
       WHERE \"Group\"='%s' AND Gene='F' AND Variation_Type='AA' AND Grouping_Type='Year'
       GROUP BY Position, AminoAcid, Year", grp))
    flagged <- scan_positions(observed_sets(d))
    data.frame(subtype = paste0("RSV:", grp), gene = "F", region = "F",
               position = as.integer(flagged), stringsAsFactors = FALSE)
  }))
}

build_covid_s <- function() {
  con <- dbConnect(duckdb::duckdb(), "data/cache/COVID/covid_explorer.duckdb", read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  # GroupType='Year_month' is one aggregation (no multi-count) and the time bucket.
  d <- dbGetQuery(con,
    "SELECT PositionBase pos, AA aa, GroupValue bucket, SUM(Count) n FROM usage
     WHERE Protein='S' AND GroupType='Year_month'
     GROUP BY PositionBase, AA, GroupValue")
  flagged <- scan_positions(observed_sets(d))
  data.frame(subtype = "COVID:SARS-CoV-2", gene = "S", region = "S",
             position = as.integer(flagged), stringsAsFactors = FALSE)
}

# ---- write -----------------------------------------------------------------

# Append glyco rows (one per subtype+region) to a file, replacing any existing
# ones. `s` holds the glyco sites for the subtypes that belong to this file.
write_file <- function(path, s) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  is_data <- !grepl("^\\s*#", lines) & nzchar(trimws(lines))
  drop <- vapply(seq_along(lines), function(i) {
    if (!is_data[i]) return(FALSE)
    f <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    length(f) >= GROUP_COL && identical(trimws(f[GROUP_COL]), GROUP)
  }, logical(1))
  lines <- lines[!drop]
  new_rows <- character(0)
  for (st in unique(s$subtype)) {
    ss <- s[s$subtype == st, ]
    for (rg in sort(unique(ss$region))) {
      r <- ss[ss$region == rg, ]
      pos <- sort(unique(r$position))
      new_rows <- c(new_rows, sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
        st, r$gene[1], rg, SITE, GROUP, paste(pos, collapse = ","), COLOR, SOURCE))
      cat(sprintf("   %-18s %-4s %2d sites: %s\n", st, rg, length(pos),
                  paste0("N", pos, collapse = ",")))
    }
  }
  writeLines(c(lines, new_rows), path)
}

sites <- rbind(build_flu_ha(), build_rsv_f(), build_covid_s())
sites <- sites[!is.na(sites$position), ]
for (path in unique(FILE_FOR)) {
  cat(basename(path), "\n")
  write_file(path, sites[FILE_FOR[sites$subtype] == path, ])
}
