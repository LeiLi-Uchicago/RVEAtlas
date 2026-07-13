#!/usr/bin/env Rscript
# =============================================================================
# generate_gene_regions.R
#
# Writes www/structures/gene_regions.tsv -- structural/functional region bands
# drawn on the Conservation (entropy) plot. Columns:
#   subtype  gene  region  start  end  color
# `start`/`end` are in the app Position coordinate (= the entropy x-axis, which
# equals the reference numbering: HA app position, Wuhan Spike numbering, RSV F0).
#
# HA regions are DERIVED from ha_numbering_review_table.csv per subtype (signal
# peptide = the run before HA1; HA1/HA2 from the table's HA_Region), so they stay
# correct per subtype. The C-terminal TM/cytoplasmic tail is the stretch past the
# HA2 ectodomain (from the table's HA2 end to the last modelled app position).
# Spike and RSV F regions are literature constants in their native numbering.
#
#   Rscript tools/generate_gene_regions.R
# =============================================================================
suppressWarnings(suppressMessages(library(DBI)))

OUT <- "www/structures/gene_regions.tsv"
NUMBERING <- "ha_numbering_review_table.csv"

# Region -> colour (rendered translucent by the plot).
COL <- c(
  "Signal peptide" = "#7f8c8d", "HA1" = "#2980b9", "HA2" = "#e67e22",
  "TM / cytoplasmic" = "#34495e", "NTD" = "#2980b9", "RBD" = "#16a085",
  "Fusion peptide" = "#c0392b", "HR1" = "#8e44ad", "HR2" = "#9b59b6",
  "TM" = "#34495e", "F2" = "#2980b9", "p27" = "#f1c40f", "F1" = "#e67e22"
)
col_of <- function(region) unname(ifelse(region %in% names(COL), COL[region], "#95a5a6"))

ha_max_position <- function() {
  db <- "data/cache/FLU/flu_explorer.duckdb"
  fallback <- c(H1N1 = 566, H3N2 = 566, H5NX = 569)
  if (!file.exists(db)) return(fallback)
  out <- tryCatch({
    con <- dbConnect(duckdb::duckdb(), db, read_only = TRUE)
    on.exit(dbDisconnect(con, shutdown = TRUE))
    vapply(names(fallback), function(st) {
      m <- dbGetQuery(con, sprintf(
        "SELECT MAX(Position) mx FROM usage WHERE \"Group\"='%s' AND Gene='HA' AND Variation_Type='AA'", st))$mx
      if (length(m) == 0 || is.na(m)) fallback[[st]] else as.integer(m)
    }, integer(1))
  }, error = function(e) fallback)
  out
}

ha_regions <- function() {
  m <- read.csv(NUMBERING, stringsAsFactors = FALSE, check.names = FALSE)
  m <- m[!as.logical(m$Is_Alignment_Gap), ]
  maxpos <- ha_max_position()
  do.call(rbind, lapply(c("H1N1", "H3N2", "H5NX"), function(st) {
    s <- m[m$Subtype == st, ]
    ha1 <- range(s$Full_HA_Position[s$HA_Region == "HA1"])
    ha2 <- range(s$Full_HA_Position[s$HA_Region == "HA2"])
    reg <- rbind(
      data.frame(region = "Signal peptide", start = 1L,          end = ha1[1] - 1L),
      data.frame(region = "HA1",            start = ha1[1],      end = ha1[2]),
      data.frame(region = "HA2",            start = ha2[1],      end = ha2[2]),
      data.frame(region = "TM / cytoplasmic", start = ha2[2] + 1L, end = maxpos[[st]])
    )
    data.frame(subtype = st, gene = "HA", reg, stringsAsFactors = FALSE)
  }))
}

# Literature constants (native = app numbering).
spike_regions <- function() {
  reg <- data.frame(
    region = c("Signal peptide", "NTD", "RBD", "Fusion peptide", "HR1", "HR2", "TM"),
    start  = c(1, 14, 319, 816, 912, 1163, 1214),
    end    = c(13, 305, 541, 837, 984, 1213, 1237), stringsAsFactors = FALSE)
  data.frame(subtype = "COVID:SARS-CoV-2", gene = "S", reg, stringsAsFactors = FALSE)
}

rsv_f_regions <- function() {
  reg <- data.frame(
    region = c("Signal peptide", "F2", "p27", "F1", "TM"),
    start  = c(1, 26, 110, 137, 525),
    end    = c(25, 109, 136, 524, 550), stringsAsFactors = FALSE)
  do.call(rbind, lapply(c("RSV:A", "RSV:B"), function(st)
    data.frame(subtype = st, gene = "F", reg, stringsAsFactors = FALSE)))
}

all_regions <- rbind(ha_regions(), spike_regions(), rsv_f_regions())
all_regions$color <- col_of(all_regions$region)
all_regions <- all_regions[, c("subtype", "gene", "region", "start", "end", "color")]

header <- c(
  "# Structural / functional regions drawn as translucent bands on the Conservation",
  "# (entropy) plot. Positions are in the app Position coordinate (= entropy x-axis).",
  "# HA rows are derived from ha_numbering_review_table.csv; Spike/RSV-F are literature",
  "# constants in native numbering. Regenerate with tools/generate_gene_regions.R.",
  paste(c("subtype", "gene", "region", "start", "end", "color"), collapse = "\t")
)
lines <- apply(all_regions, 1, function(r) paste(trimws(r), collapse = "\t"))
writeLines(c(header, lines), OUT)
cat("Wrote", OUT, "with", nrow(all_regions), "region rows:\n")
print(all_regions, row.names = FALSE)
