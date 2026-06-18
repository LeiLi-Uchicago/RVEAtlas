conservation_effective_group <- function(filter_enabled, selected_group) {
  if (!isTRUE(filter_enabled)) return("All")
  if (is.null(selected_group) || length(selected_group) == 0) return("All")
  selected_group <- as.character(selected_group[[1]])
  if (!nzchar(selected_group)) "All" else selected_group
}
