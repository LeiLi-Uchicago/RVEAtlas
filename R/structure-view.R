# =============================================================================
# structure-view.R
#
# Helpers for the 3D protein-structure views (r3dmol / 3Dmol.js) used on the
# Conservation (entropy) and Single Site pages.
#
# Design: the transform functions are pure (no Shiny, no global state) so they
# can be unit-tested; file loaders are thin wrappers around them. Widget /
# proxy calls live in server.R.
#
# Position model (H1N1 HA, PDB 4JTV; see memory 4jtv-ha-numbering-mapping):
#   app Full_HA_Position --(HA numbering table)--> HA_Region + Numbering_Position
#   HA_Region -> PDB chain (HA1 -> A, HA2 -> B); resi = Numbering_Position.
# Other pathogens/genes slot in by adding rows to structure_config.tsv.
# =============================================================================

# ---- Locations -------------------------------------------------------------

# read.delim with comment.char = "#" would also strip in-field values that
# start with "#" (e.g. hex colors like #4363d8). Instead drop only full-line
# comments (first non-space char is "#") and parse the rest.
sv_read_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  keep <- !grepl("^\\s*#", lines)
  lines <- lines[keep]
  if (length(lines) == 0) return(data.frame())
  # na.strings = "" (not the default "NA") so the literal gene/region value "NA"
  # (neuraminidase) is kept as a string rather than read as a missing value.
  utils::read.delim(
    text = paste(lines, collapse = "\n"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = ""
  )
}

sv_structure_dir <- function() {
  cand <- c(
    Sys.getenv("RVEATLAS_STRUCTURE_DIR", ""),
    file.path("www", "structures")
  )
  cand <- cand[nzchar(cand)]
  hit <- cand[dir.exists(cand)]
  if (length(hit) > 0) hit[[1]] else file.path("www", "structures")
}

# ---- Structure config ------------------------------------------------------

sv_empty_config <- function() {
  data.frame(
    pathogen = character(), subtype = character(), gene = character(),
    pdb_id = character(), structure_file = character(),
    region_chains = character(), resi_from_numbering = logical(),
    epitope_file = character(), title = character(),
    stringsAsFactors = FALSE
  )
}

sv_load_structure_config <- function(dir = sv_structure_dir()) {
  path <- file.path(dir, "structure_config.tsv")
  if (!file.exists(path)) return(sv_empty_config())
  df <- sv_read_tsv(path)
  if (!is.null(df$resi_from_numbering)) {
    df$resi_from_numbering <- as.logical(df$resi_from_numbering)
  }
  df
}

# A subtype/gene can have several structures (config rows); `pdb_id` picks one,
# otherwise the first row is the default.
sv_get_structure_config <- function(subtype, gene, config = sv_load_structure_config(),
                                     pdb_id = NULL) {
  if (is.null(subtype) || is.null(gene) || nrow(config) == 0) return(NULL)
  hit <- config[
    config$subtype == as.character(subtype) & config$gene == as.character(gene),
    , drop = FALSE
  ]
  if (nrow(hit) == 0) return(NULL)
  if (!is.null(pdb_id) && length(pdb_id) == 1 && !is.na(pdb_id) && nzchar(pdb_id)) {
    pick <- hit[as.character(hit$pdb_id) == as.character(pdb_id), , drop = FALSE]
    if (nrow(pick) > 0) hit <- pick
  }
  as.list(hit[1, , drop = FALSE])
}

# All structures available for a subtype/gene, as data.frame(pdb_id, label) in
# config order (first = default). `label` falls back to pdb_id when absent.
sv_list_structures <- function(subtype, gene, config = sv_load_structure_config()) {
  empty <- data.frame(pdb_id = character(0), label = character(0), stringsAsFactors = FALSE)
  if (is.null(subtype) || is.null(gene) || nrow(config) == 0) return(empty)
  hit <- config[
    config$subtype == as.character(subtype) & config$gene == as.character(gene),
    , drop = FALSE
  ]
  if (nrow(hit) == 0) return(empty)
  pdb <- as.character(hit$pdb_id)
  lab <- if (!is.null(hit$label)) as.character(hit$label) else rep(NA_character_, length(pdb))
  blank <- is.na(lab) | !nzchar(lab)
  lab[blank] <- pdb[blank]
  data.frame(pdb_id = pdb, label = lab, stringsAsFactors = FALSE)
}

# "HA1:A;HA2:B"       -> list(HA1 = "A", HA2 = "B")
# "HA1:A,C,E;HA2:B,D,F"-> list(HA1 = c("A","C","E"), HA2 = c("B","D","F"))
# One region can map to several chains (e.g. the three protomers of an HA
# trimer), so each region resolves to a vector of chain IDs.
sv_parse_region_chains <- function(str) {
  if (is.null(str) || length(str) == 0 || is.na(str) || !nzchar(str)) {
    return(stats::setNames(list(), character(0)))
  }
  parts <- trimws(strsplit(str, ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  kv <- strsplit(parts, ":", fixed = TRUE)
  regions <- vapply(kv, function(x) trimws(x[[1]]), character(1))
  chains <- lapply(kv, function(x) {
    if (length(x) < 2) return(character(0))
    ch <- trimws(strsplit(x[[2]], ",", fixed = TRUE)[[1]])
    ch[nzchar(ch)]
  })
  stats::setNames(chains, regions)
}

# ---- HA numbering ----------------------------------------------------------

sv_empty_numbering <- function() {
  data.frame(
    Subtype = character(), Full_HA_Position = numeric(),
    HA_Region = character(), Numbering_Position = numeric(),
    stringsAsFactors = FALSE
  )
}

# Identity numbering for structures whose PDB residues already equal the app
# position (resi_from_numbering = FALSE), e.g. NA tetramers renumbered by
# tools/renumber_na.py. Every app position maps to itself under a single region
# (the config's region_chains region, e.g. "NA"). Used in place of the HA
# numbering table so genes without an alignment table still highlight residues.
sv_identity_numbering <- function(region, positions) {
  region <- region[!is.na(region) & nzchar(region)]
  positions <- suppressWarnings(as.numeric(positions))
  positions <- sort(unique(positions[!is.na(positions)]))
  if (length(region) == 0 || length(positions) == 0) return(sv_empty_numbering())
  do.call(rbind, lapply(region, function(rg) {
    data.frame(
      Subtype = NA_character_,
      Full_HA_Position = positions,
      HA_Region = rg,
      Numbering_Position = positions,
      stringsAsFactors = FALSE
    )
  }))
}

# Residue numbers actually present in a structure's protein chains. Used to build
# identity numbering so that positions not modeled in the crystal (e.g. the NA
# stalk/transmembrane region, or disordered loops) are correctly reported as
# "not resolved" rather than silently mapped. Sentinel numbers >= 1000 written by
# the renumber tools for unalignable structure residues are excluded.
sv_structure_resolved_positions <- function(structure_file, dir = sv_structure_dir()) {
  path <- file.path(dir, structure_file)
  if (!file.exists(path)) return(numeric(0))
  lines <- readLines(path, warn = FALSE)
  atom <- lines[substr(lines, 1, 6) == "ATOM  "]
  if (length(atom) == 0) return(numeric(0))
  resseq <- suppressWarnings(as.numeric(substr(atom, 23, 26)))
  resseq <- resseq[!is.na(resseq) & resseq < 1000]
  sort(unique(resseq))
}

sv_load_numbering <- function(subtype = NULL, path = NULL) {
  if (is.null(path)) {
    path <- if (exists("ha_numbering_source_path", mode = "function")) {
      ha_numbering_source_path()
    } else {
      "ha_numbering_review_table.csv"
    }
  }
  if (length(path) == 0 || is.na(path) || !file.exists(path)) return(sv_empty_numbering())
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("Subtype", "Full_HA_Position", "HA_Region", "Numbering_Position", "Is_Alignment_Gap")
  if (!all(req %in% names(m))) return(sv_empty_numbering())
  is_gap <- as.logical(m$Is_Alignment_Gap)
  is_gap[is.na(is_gap)] <- FALSE
  m <- m[!is_gap & !is.na(m$Full_HA_Position), , drop = FALSE]
  out <- data.frame(
    Subtype = as.character(m$Subtype),
    Full_HA_Position = suppressWarnings(as.numeric(m$Full_HA_Position)),
    HA_Region = as.character(m$HA_Region),
    Numbering_Position = suppressWarnings(as.numeric(m$Numbering_Position)),
    stringsAsFactors = FALSE
  )
  if (!is.null(subtype)) out <- out[out$Subtype == as.character(subtype), , drop = FALSE]
  out
}

# ---- Position mapping (pure) ----------------------------------------------

# app Full_HA_Position vector -> data.frame(position, chain, resi) for the
# residues that map cleanly; unmapped positions are dropped.
sv_map_positions <- function(positions, numbering_df, region_chains) {
  positions <- suppressWarnings(as.numeric(positions))
  empty <- data.frame(position = numeric(), chain = character(), resi = numeric(),
                      stringsAsFactors = FALSE)
  if (length(positions) == 0 || nrow(numbering_df) == 0 || length(region_chains) == 0) {
    return(empty)
  }
  idx <- match(positions, numbering_df$Full_HA_Position)
  region <- numbering_df$HA_Region[idx]
  resi <- numbering_df$Numbering_Position[idx]
  keep <- !is.na(idx) & !is.na(resi) & !is.na(region) & region %in% names(region_chains)
  if (!any(keep)) return(empty)
  pos_k <- positions[keep]; reg_k <- region[keep]; resi_k <- resi[keep]
  rows <- lapply(seq_along(pos_k), function(i) {
    chains <- region_chains[[reg_k[i]]]
    chains <- chains[!is.na(chains) & nzchar(chains)]
    if (length(chains) == 0) return(NULL)
    data.frame(position = pos_k[i], chain = chains, resi = resi_k[i],
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

# ---- Entropy -> colors (pure) ---------------------------------------------

# Low entropy (conserved) -> first palette color; high entropy (variable) ->
# last. Returns per-value color/bin plus a legend table.
sv_entropy_colors <- function(entropy, n_bins = 9, domain = NULL,
                              palette = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) {
  entropy <- suppressWarnings(as.numeric(entropy))
  if (is.null(domain)) {
    dmax <- suppressWarnings(max(entropy, na.rm = TRUE))
    if (!is.finite(dmax) || dmax <= 0) dmax <- 1
    domain <- c(0, dmax)
  }
  n_bins <- max(1L, as.integer(n_bins))
  breaks <- seq(domain[1], domain[2], length.out = n_bins + 1)
  cols <- grDevices::colorRampPalette(palette)(n_bins)
  bin <- findInterval(entropy, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  bin[is.na(entropy)] <- NA_integer_
  color <- ifelse(is.na(bin), NA_character_, cols[bin])
  legend <- data.frame(
    bin = seq_len(n_bins),
    lower = breaks[-(n_bins + 1)],
    upper = breaks[-1],
    color = cols,
    stringsAsFactors = FALSE
  )
  list(color = color, bin = bin, breaks = breaks, colors = cols,
       legend = legend, domain = domain)
}

# Combine an entropy table (columns Position, Entropy) with the numbering map
# into per-residue colors ready for the viewer.
sv_entropy_residue_colors <- function(entropy_df, numbering_df, region_chains,
                                       n_bins = 9, domain = NULL,
                                       palette = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) {
  empty <- data.frame(position = numeric(), chain = character(), resi = numeric(),
                      entropy = numeric(), color = character(), bin = integer(),
                      stringsAsFactors = FALSE)
  if (is.null(entropy_df) || nrow(entropy_df) == 0) {
    return(list(residues = empty, legend = sv_entropy_colors(numeric(0))$legend, domain = c(0, 1)))
  }
  mapped <- sv_map_positions(entropy_df$Position, numbering_df, region_chains)
  if (nrow(mapped) == 0) {
    return(list(residues = empty, legend = sv_entropy_colors(entropy_df$Entropy, n_bins = n_bins, domain = domain, palette = palette)$legend, domain = c(0, 1)))
  }
  ent <- entropy_df$Entropy[match(mapped$position, entropy_df$Position)]
  col <- sv_entropy_colors(ent, n_bins = n_bins, domain = domain, palette = palette)
  mapped$entropy <- ent
  mapped$color <- col$color
  mapped$bin <- col$bin
  mapped <- mapped[!is.na(mapped$color), , drop = FALSE]
  list(residues = mapped, legend = col$legend, domain = col$domain)
}

# ---- Epitopes --------------------------------------------------------------

sv_empty_epitopes <- function() {
  data.frame(site = character(), region = character(), position = numeric(),
             color = character(), group = character(), stringsAsFactors = FALSE)
}

# Loads the epitope table referenced by a structure-config row and expands the
# comma-separated `positions` field into one row per residue. Preserves `group`
# column if present for filtering by epitope group.
sv_load_epitopes <- function(config_row, subtype, gene, dir = sv_structure_dir()) {
  if (is.null(config_row)) return(sv_empty_epitopes())
  ef <- config_row$epitope_file
  if (is.null(ef) || length(ef) == 0 || is.na(ef) || !nzchar(ef)) return(sv_empty_epitopes())
  path <- file.path(dir, ef)
  if (!file.exists(path)) return(sv_empty_epitopes())
  df <- sv_read_tsv(path)
  need <- c("subtype", "gene", "region", "site", "positions", "color")
  if (!all(need %in% names(df))) return(sv_empty_epitopes())
  has_group <- "group" %in% names(df)
  df <- df[df$subtype == as.character(subtype) & df$gene == as.character(gene), , drop = FALSE]
  if (nrow(df) == 0) return(sv_empty_epitopes())
  rows <- lapply(seq_len(nrow(df)), function(i) {
    pos <- suppressWarnings(as.numeric(strsplit(gsub("[[:space:]]", "", df$positions[i]), ",", fixed = TRUE)[[1]]))
    pos <- pos[!is.na(pos)]
    if (length(pos) == 0) return(NULL)
    res <- data.frame(site = df$site[i], region = df$region[i], position = pos,
                      color = df$color[i], stringsAsFactors = FALSE)
    if (has_group) res$group <- df$group[i]
    res
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(sv_empty_epitopes())
  do.call(rbind, rows)
}

# List all unique epitope groups in a loaded epitope table.
sv_epitope_groups <- function(epitope_df) {
  if (is.null(epitope_df) || nrow(epitope_df) == 0 || !"group" %in% names(epitope_df)) {
    return(character(0))
  }
  sort(unique(epitope_df$group[!is.na(epitope_df$group) & nzchar(epitope_df$group)]))
}

# Filter epitope table by group(s). If groups is NULL/empty, all epitopes returned.
sv_filter_epitopes_by_group <- function(epitope_df, groups = NULL) {
  if (is.null(epitope_df) || nrow(epitope_df) == 0) return(epitope_df)
  if (is.null(groups) || length(groups) == 0) return(epitope_df)
  if (!"group" %in% names(epitope_df)) return(epitope_df)
  epitope_df[epitope_df$group %in% groups, , drop = FALSE]
}

# Epitope table -> viewer residues. Epitope positions are already in numbering
# space (= PDB resi when resi_from_numbering is TRUE); region -> chain.
# Preserves `group` column if present.
sv_epitope_residues <- function(epitope_df, region_chains) {
  empty <- data.frame(site = character(), chain = character(), resi = numeric(),
                      color = character(), group = character(), stringsAsFactors = FALSE)
  if (is.null(epitope_df) || nrow(epitope_df) == 0 || length(region_chains) == 0) return(empty)
  keep <- !is.na(epitope_df$position) & epitope_df$region %in% names(region_chains)
  if (!any(keep)) return(empty)
  ep <- epitope_df[keep, , drop = FALSE]
  has_group <- "group" %in% names(ep)
  rows <- lapply(seq_len(nrow(ep)), function(i) {
    chains <- region_chains[[ep$region[i]]]
    chains <- chains[!is.na(chains) & nzchar(chains)]
    if (length(chains) == 0) return(NULL)
    res <- data.frame(site = ep$site[i], chain = chains, resi = ep$position[i],
                      color = ep$color[i], stringsAsFactors = FALSE)
    if (has_group) res$group <- ep$group[i]
    res
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

# Distinct (site, color) legend rows, preserving first-seen order.
sv_epitope_legend <- function(epitope_df) {
  if (is.null(epitope_df) || nrow(epitope_df) == 0) {
    return(data.frame(site = character(), color = character(), stringsAsFactors = FALSE))
  }
  keys <- !duplicated(epitope_df$site)
  data.frame(site = epitope_df$site[keys], color = epitope_df$color[keys],
             stringsAsFactors = FALSE)
}

# Find which epitope group(s) contain a specific position in a given region.
# Returns character vector of group names, or empty if position is not in any epitope.
sv_position_epitope_groups <- function(position, region, epitope_df) {
  if (is.null(epitope_df) || nrow(epitope_df) == 0) return(character(0))
  if (!"group" %in% names(epitope_df)) return(character(0))
  if (is.na(position) || is.na(region)) return(character(0))
  hits <- epitope_df[epitope_df$position == position & epitope_df$region == region, ]
  if (nrow(hits) == 0) return(character(0))
  sort(unique(hits$group[!is.na(hits$group) & nzchar(hits$group)]))
}

# ---- Viewer builders (r3dmol) ----------------------------------------------
# These reference r3dmol only inside the function bodies, so the module still
# sources cleanly (for unit tests) without the package attached.

sv_pdb_text <- function(config_row, dir = sv_structure_dir()) {
  if (is.null(config_row)) return(NULL)
  path <- file.path(dir, config_row$structure_file)
  if (!file.exists(path)) return(NULL)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

sv_view_modes <- function() {
  c("Surface" = "surface", "Cartoon" = "cartoon", "Stick" = "stick", "Sphere" = "sphere")
}

# Backbone/atom representation for a non-surface mode.
sv_style_for_mode <- function(mode, color) {
  switch(mode,
    stick  = r3dmol::m_style_stick(color = color),
    sphere = r3dmol::m_style_sphere(color = color),
    r3dmol::m_style_cartoon(color = color)
  )
}

# Solid-color surface style. m_style_surface() injects colorscheme = "default"
# (CPK element coloring), which overrides `color` in 3Dmol; building the spec
# without colorscheme makes the surface honor the explicit `color`.
sv_surface_style <- function(opacity = 0.9, color = NULL) {
  spec <- list(opacity = opacity)
  if (!is.null(color)) spec$color <- color
  structure(spec, class = "SurfaceStyleSpec")
}

# Unique (chain, resi) residues in a PDB (one CA per residue).
sv_pdb_residues <- function(pdb_text) {
  lines <- strsplit(pdb_text, "\n", fixed = TRUE)[[1]]
  atom <- lines[substr(lines, 1, 4) == "ATOM"]
  ca <- atom[substr(atom, 13, 16) == " CA "]
  if (length(ca) == 0) return(data.frame(chain = character(), resi = integer(), stringsAsFactors = FALSE))
  df <- data.frame(
    chain = substr(ca, 22, 22),
    resi = suppressWarnings(as.integer(substr(ca, 23, 26))),
    stringsAsFactors = FALSE
  )
  df[!is.na(df$resi), , drop = FALSE]
}

# Assign every residue exactly one color: base_color, then `residues` overrides,
# then the selected residue. Disjoint groups avoid overlapping surfaces (which
# blend into muddy colors) and the CPK default.
sv_residue_color_map <- function(pdb_text, residues, base_color, current, current_color) {
  allres <- sv_pdb_residues(pdb_text)
  if (nrow(allres) == 0) return(allres)
  allres$color <- base_color
  keyall <- paste(allres$chain, allres$resi)
  if (!is.null(residues) && nrow(residues) > 0) {
    idx <- match(paste(residues$chain, residues$resi), keyall)
    ok <- !is.na(idx)
    allres$color[idx[ok]] <- residues$color[ok]
  }
  if (!is.null(current) && nrow(current) > 0) {
    idx <- match(paste(current$chain, current$resi), keyall)
    ok <- !is.na(idx)
    allres$color[idx[ok]] <- current_color
  }
  allres
}

# Unified viewer builder shared by both pages.
#   residues       data.frame(chain, resi, color) to color (entropy bins or epitopes)
#   mode           "surface" | "cartoon" | "stick" | "sphere"
#   current        optional data.frame(chain, resi) for the highlighted residue
#   current_label  optional text drawn as a label on the highlighted residue
# In surface mode every residue gets a disjoint, explicitly colored surface
# patch; in the other modes the same colors are applied to that representation.
# The highlight therefore uses the same representation as the rest of the view.
sv_build_viewer <- function(pdb_text, residues = NULL, mode = "surface",
                            base_color = "#c8ced6", background = "#ffffff",
                            current = NULL, current_color = "#ff2d2d",
                            current_label = NULL, surface_opacity = 0.92) {
  is_surface <- identical(mode, "surface")

  v <- r3dmol::r3dmol(backgroundColor = background) %>%
    r3dmol::m_add_model(data = pdb_text, format = "pdb")

  if (is_surface) {
    # Thin cartoon underneath for structural context, then colored surface patches.
    v <- v %>% r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = base_color))
    allres <- sv_residue_color_map(pdb_text, residues, base_color, current, current_color)
    key <- paste(allres$chain, allres$color, sep = "|")
    for (k in unique(key)) {
      sub <- allres[key == k, , drop = FALSE]
      v <- v %>% r3dmol::m_add_surface(
        type = "VDW",
        style = sv_surface_style(surface_opacity, color = sub$color[[1]]),
        atomsel = r3dmol::m_sel(chain = sub$chain[[1]], resi = sub$resi)
      )
    }
  } else {
    v <- v %>% r3dmol::m_set_style(style = sv_style_for_mode(mode, base_color))
    if (!is.null(residues) && nrow(residues) > 0) {
      key <- paste(residues$chain, residues$color, sep = "|")
      for (k in unique(key)) {
        sub <- residues[key == k, , drop = FALSE]
        v <- v %>% r3dmol::m_set_style(
          sel = r3dmol::m_sel(chain = sub$chain[[1]], resi = sub$resi),
          style = sv_style_for_mode(mode, sub$color[[1]])
        )
      }
    }
    if (!is.null(current) && nrow(current) > 0) {
      for (ch in unique(current$chain)) {
        sub <- current[current$chain == ch, , drop = FALSE]
        v <- v %>% r3dmol::m_set_style(
          sel = r3dmol::m_sel(chain = ch, resi = sub$resi),
          style = sv_style_for_mode(mode, current_color)
        )
      }
    }
  }

  if (!is.null(current) && nrow(current) > 0 &&
      !is.null(current_label) && nzchar(current_label)) {
    # One label only, on the first protomer, to avoid three overlapping tags.
    label_chain <- current$chain[[1]]
    label_resi <- current$resi[current$chain == label_chain]
    v <- v %>% r3dmol::m_add_label(
      text = current_label,
      sel = r3dmol::m_sel(chain = label_chain, resi = label_resi),
      style = r3dmol::m_style_label(
        font = "Arial", fontSize = 18, fontColor = "white",
        borderThickness = 1.2, borderColor = "#1a1a1a", borderOpacity = 1,
        backgroundColor = current_color, backgroundOpacity = 0.95,
        inFront = TRUE
      )
    )
  }

  v %>% r3dmol::m_zoom_to()
}
