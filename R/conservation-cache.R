CONSERVATION_CACHE_SCHEMA_VERSION <- 1L
CONSERVATION_CACHE_REQUIRED_COLUMNS <- c(
  "Subtype", "Variation_Type", "Gene", "Position", "Entropy", "Pos_Total"
)

conservation_entropy_from_counts <- function(counts) {
  required <- c("Position", "AminoAcid", "Count")
  if (is.null(counts) || !all(required %in% names(counts)) || nrow(counts) == 0) {
    return(data.frame(
      Position = numeric(),
      Entropy = numeric(),
      Pos_Total = numeric()
    ))
  }

  counts <- counts[
    !is.na(counts$AminoAcid) &
      !counts$AminoAcid %in% c("X", "-") &
      !is.na(counts$Count),
    required,
    drop = FALSE
  ]
  if (nrow(counts) == 0) {
    return(data.frame(
      Position = numeric(),
      Entropy = numeric(),
      Pos_Total = numeric()
    ))
  }

  residue_counts <- stats::aggregate(
    counts$Count,
    by = list(Position = counts$Position, AminoAcid = counts$AminoAcid),
    FUN = sum,
    na.rm = TRUE
  )
  names(residue_counts)[[3]] <- "AA_Sum"

  split_counts <- split(residue_counts, residue_counts$Position, drop = TRUE)
  rows <- lapply(split_counts, function(position_counts) {
    total <- sum(position_counts$AA_Sum, na.rm = TRUE)
    probabilities <- position_counts$AA_Sum / total
    probabilities <- probabilities[is.finite(probabilities) & probabilities > 0]
    data.frame(
      Position = position_counts$Position[[1]],
      Entropy = if (length(probabilities) == 0) 0 else -sum(probabilities * log2(probabilities)),
      Pos_Total = total
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(suppressWarnings(as.numeric(as.character(result$Position))), result$Position), , drop = FALSE]
}

conservation_cache_structure_valid <- function(
  cache,
  schema_version = CONSERVATION_CACHE_SCHEMA_VERSION
) {
  is.list(cache) &&
    identical(cache$schema_version, as.integer(schema_version)) &&
    is.data.frame(cache$data) &&
    all(CONSERVATION_CACHE_REQUIRED_COLUMNS %in% names(cache$data))
}

conservation_signature_identity <- function(signature) {
  if (!is.list(signature)) return(NULL)

  files <- as.character(signature$files)
  sizes <- suppressWarnings(as.numeric(signature$sizes))
  if (length(files) == 0 || length(files) != length(sizes)) return(NULL)

  list(
    files = gsub("\\\\", "/", files),
    sizes = sizes
  )
}

conservation_cache_valid <- function(
  cache,
  source_signature,
  schema_version = CONSERVATION_CACHE_SCHEMA_VERSION
) {
  if (!conservation_cache_structure_valid(cache, schema_version)) return(FALSE)

  # A deployed copy may intentionally contain only the pre-calculated cache.
  # In that case there is no source data against which to mark it stale.
  if (length(source_signature$files) == 0) return(TRUE)

  # File mtimes commonly change when the application is copied to another
  # machine. Compare the portable parts of the signature instead.
  identical(
    conservation_signature_identity(cache$source_signature),
    conservation_signature_identity(source_signature)
  )
}

conservation_cache_lookup <- function(cache, subtype, variation_type, gene) {
  if (!is.list(cache) || !is.data.frame(cache$data)) return(NULL)
  data <- cache$data
  if (!all(CONSERVATION_CACHE_REQUIRED_COLUMNS %in% names(data))) return(NULL)

  result <- data[
    as.character(data$Subtype) == as.character(subtype) &
      as.character(data$Variation_Type) == as.character(variation_type) &
      as.character(data$Gene) == as.character(gene),
    c("Position", "Entropy", "Pos_Total"),
    drop = FALSE
  ]
  if (nrow(result) == 0) NULL else result
}

conservation_cache_path <- function(pathogen_id) {
  file.path("data", "cache", pathogen_id, "conservation_entropy.rds")
}

conservation_source_files <- function(pathogen_id) {
  if (identical(pathogen_id, "FLU")) {
    candidates <- c(RDS_CACHE, DUCKDB_CACHE, DUCKDB_META_CACHE)
  } else {
    cfg <- PATHOGEN_ADAPTERS[[pathogen_id]]
    candidates <- c(cfg$duckdb, cfg$metadata)
  }
  unique(candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)])
}

conservation_source_signature <- function(pathogen_id) {
  files <- sort(conservation_source_files(pathogen_id))
  info <- file.info(files)
  list(
    files = files,
    sizes = unname(as.numeric(info$size)),
    mtimes = unname(as.numeric(info$mtime))
  )
}

conservation_cache_env <- new.env(parent = emptyenv())

build_conservation_entropy_cache <- function(pathogen_id) {
  cfg <- PATHOGEN_ADAPTERS[[pathogen_id]]
  if (is.null(cfg)) stop("Unsupported pathogen: ", pathogen_id, call. = FALSE)

  subtype_values <- unname(cfg$subtype_choices)
  result_rows <- list()
  row_index <- 0L

  for (subtype in subtype_values) {
    for (variation_type in c("AA", "NT")) {
      genes <- usage_available_genes(subtype, variation_type)
      if (length(genes) == 0) next

      for (gene in genes) {
        groups <- usage_available_groups(subtype, variation_type, gene)
        if (length(groups) == 0) next

        entropy <- tryCatch(
          usage_entropy_data(subtype, variation_type, gene, groups[[1]], "All"),
          error = function(error) {
            warning(
              "Could not precompute Conservation entropy for ",
              pathogen_id, " / ", subtype, " / ", variation_type, " / ", gene,
              ": ", conditionMessage(error),
              call. = FALSE
            )
            NULL
          }
        )
        if (is.null(entropy) || nrow(entropy) == 0) next

        row_index <- row_index + 1L
        result_rows[[row_index]] <- data.frame(
          Subtype = as.character(subtype),
          Variation_Type = as.character(variation_type),
          Gene = as.character(gene),
          Position = entropy$Position,
          Entropy = as.numeric(entropy$Entropy),
          Pos_Total = as.numeric(entropy$Pos_Total),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  data <- if (length(result_rows) == 0) {
    data.frame(
      Subtype = character(),
      Variation_Type = character(),
      Gene = character(),
      Position = numeric(),
      Entropy = numeric(),
      Pos_Total = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, result_rows)
  }

  cache <- list(
    schema_version = CONSERVATION_CACHE_SCHEMA_VERSION,
    source_signature = conservation_source_signature(pathogen_id),
    built_at = Sys.time(),
    data = data
  )
  path <- conservation_cache_path(pathogen_id)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, path)
  conservation_cache_env[[pathogen_id]] <- cache
  path
}

ensure_conservation_entropy_cache <- function(pathogen_id, force = FALSE) {
  path <- conservation_cache_path(pathogen_id)
  signature <- conservation_source_signature(pathogen_id)
  cache <- if (!isTRUE(force) && file.exists(path)) {
    tryCatch(readRDS(path), error = function(error) NULL)
  } else {
    NULL
  }

  if (!isTRUE(force) && conservation_cache_valid(cache, signature)) {
    conservation_cache_env[[pathogen_id]] <- cache
    return(path)
  }

  tryCatch(
    build_conservation_entropy_cache(pathogen_id),
    error = function(error) {
      if (!isTRUE(force) && conservation_cache_structure_valid(cache)) {
        warning(
          "Conservation entropy cache could not be refreshed for ",
          pathogen_id, "; using the existing pre-calculated cache: ",
          conditionMessage(error),
          call. = FALSE
        )
        conservation_cache_env[[pathogen_id]] <- cache
        return(path)
      }
      stop(error)
    }
  )
}

ensure_all_conservation_entropy_caches <- function(force = FALSE) {
  pathogen_ids <- available_pathogen_ids()
  invisible(lapply(pathogen_ids, function(pathogen_id) {
    tryCatch(
      ensure_conservation_entropy_cache(pathogen_id, force = force),
      error = function(error) {
        warning(
          "Conservation entropy cache could not be refreshed for ",
          pathogen_id, ": ", conditionMessage(error),
          call. = FALSE
        )
        NULL
      }
    )
  }))
}

load_conservation_entropy_cache <- function(pathogen_id) {
  cached <- conservation_cache_env[[pathogen_id]]
  if (!is.null(cached)) return(cached)

  path <- tryCatch(
    ensure_conservation_entropy_cache(pathogen_id),
    error = function(error) conservation_cache_path(pathogen_id)
  )
  if (is.null(path) || !file.exists(path)) return(NULL)
  cache <- tryCatch(readRDS(path), error = function(error) NULL)
  if (!conservation_cache_structure_valid(cache)) return(NULL)
  conservation_cache_env[[pathogen_id]] <- cache
  cache
}

conservation_cached_entropy <- function(subtype, variation_type, gene) {
  pathogen_id <- pathogen_from_subtype(subtype)
  if (is.na(pathogen_id)) return(NULL)
  conservation_cache_lookup(
    load_conservation_entropy_cache(pathogen_id),
    subtype,
    variation_type,
    gene
  )
}
