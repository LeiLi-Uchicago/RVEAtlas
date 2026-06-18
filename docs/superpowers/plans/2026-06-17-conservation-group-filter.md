# Conservation Group Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the Conservation gene selector left, place optional grouping controls in a disabled-by-default right-side box, and warn users whenever they enable the slower real-time group filter.

**Architecture:** Add one pure helper for resolving the effective group (`"All"` while disabled, selected group while enabled), then use it everywhere Conservation calculations and labels choose a group. Keep the existing Shiny inputs and data-loading flow, adding only a Boolean switch, control-state observer, modal, responsive styling, and selection-preserving dropdown updates.

**Tech Stack:** R, Shiny, shinyjs, bslib, testthat, CSS, in-app browser

---

## File Structure

- Create `R/conservation-filter.R`: pure effective-group helper, isolated for unit testing.
- Create `tests/testthat.R`: repository test entry point.
- Create `tests/testthat/test-conservation-filter.R`: behavior tests for disabled/enabled group resolution.
- Create `tests/testthat/test-conservation-filter-ui.R`: focused source-level contract checks for required UI/server wiring.
- Modify `global.R`: source the new helper during app startup.
- Modify `ui.R`: approved layout A, switch markup, disabled styling, and responsive stacking.
- Modify `server.R`: enable/disable controls, modal warning, preserve selections, and use the effective group.

### Task 1: Effective Group Behavior

**Files:**
- Create: `R/conservation-filter.R`
- Create: `tests/testthat.R`
- Create: `tests/testthat/test-conservation-filter.R`
- Modify: `global.R:1-30`

- [ ] **Step 1: Write the failing helper tests**

```r
source(file.path("R", "conservation-filter.R"), local = TRUE)

testthat::test_that("disabled Conservation filtering always uses all groups", {
  testthat::expect_identical(
    conservation_effective_group(FALSE, "3C.2a1b"),
    "All"
  )
})

testthat::test_that("enabled Conservation filtering uses the selected group", {
  testthat::expect_identical(
    conservation_effective_group(TRUE, "3C.2a1b"),
    "3C.2a1b"
  )
})

testthat::test_that("enabled filtering safely falls back when no group is selected", {
  testthat::expect_identical(conservation_effective_group(TRUE, NULL), "All")
  testthat::expect_identical(conservation_effective_group(TRUE, ""), "All")
})
```

Add this runner to `tests/testthat.R`:

```r
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The testthat package is required to run tests.")
}

testthat::test_dir("tests/testthat")
```

- [ ] **Step 2: Run the tests and verify RED**

Run:

```bash
Rscript tests/testthat.R
```

Expected: FAIL because `conservation_effective_group()` does not exist.

- [ ] **Step 3: Implement the minimal helper**

Create `R/conservation-filter.R`:

```r
conservation_effective_group <- function(filter_enabled, selected_group) {
  if (!isTRUE(filter_enabled)) return("All")
  if (is.null(selected_group) || length(selected_group) == 0) return("All")

  selected_group <- as.character(selected_group[[1]])
  if (!nzchar(selected_group)) "All" else selected_group
}
```

Source it near the top of `global.R`:

```r
source(file.path("R", "conservation-filter.R"), local = TRUE)
```

- [ ] **Step 4: Run the helper tests and verify GREEN**

Run:

```bash
Rscript tests/testthat.R
```

Expected: all three helper tests PASS.

- [ ] **Step 5: Commit the helper**

```bash
git add R/conservation-filter.R tests/testthat.R tests/testthat/test-conservation-filter.R global.R
git commit -m "test: define Conservation filter behavior"
```

### Task 2: Conservation Control Layout

**Files:**
- Create: `tests/testthat/test-conservation-filter-ui.R`
- Modify: `ui.R:470-590`
- Modify: `ui.R:1104-1120`

- [ ] **Step 1: Write the failing UI contract test**

```r
testthat::test_that("Conservation controls expose the approved filter box", {
  ui_source <- paste(readLines("ui.R", warn = FALSE), collapse = "\n")

  testthat::expect_match(
    ui_source,
    'checkboxInput\\("ent_filter_enabled",\\s*"Filter by group",\\s*value = FALSE\\)'
  )
  testthat::expect_match(ui_source, 'class = "ent-filter-box"')
  testthat::expect_match(ui_source, 'class = "ent-filter-fields"')

  gene_position <- regexpr('selectInput\\("ent_gene"', ui_source)[[1]]
  filter_position <- regexpr('class = "ent-filter-box"', ui_source)[[1]]
  testthat::expect_gt(gene_position, 0)
  testthat::expect_gt(filter_position, gene_position)
})

testthat::test_that("Conservation filter fields stack responsively", {
  ui_source <- paste(readLines("ui.R", warn = FALSE), collapse = "\n")

  testthat::expect_match(
    ui_source,
    '@media \\(max-width: 768px\\)[\\s\\S]*\\.ent-filter-fields[\\s\\S]*grid-template-columns:\\s*1fr'
  )
})
```

- [ ] **Step 2: Run the tests and verify RED**

Run:

```bash
Rscript tests/testthat.R
```

Expected: the UI contract tests FAIL because the switch, box, and responsive layout do not exist.

- [ ] **Step 3: Implement layout A and switch styling**

Replace the first Conservation control row with:

```r
fluidRow(
  column(
    3,
    helpText("Calculates Shannon Entropy to identify highly conserved valleys and hypervariable peaks across the entire gene. Subtype is controlled globally.")
  ),
  column(
    3,
    div(
      class = "ent-gene-control",
      selectInput("ent_gene", "Gene:", choices = NULL)
    )
  ),
  column(
    6,
    div(
      class = "ent-filter-box ent-filter-disabled",
      div(
        class = "ent-filter-header",
        checkboxInput("ent_filter_enabled", "Filter by group", value = FALSE)
      ),
      div(
        class = "ent-filter-fields",
        selectInput("ent_group_by", "Group by:", choices = NULL),
        selectInput("ent_group", "Group:", choices = NULL)
      )
    )
  )
)
```

Add focused CSS to the existing `tags$style()` block:

```css
.ent-filter-box {
  border: 1px solid var(--rve-border);
  border-radius: var(--rve-radius);
  background: var(--rve-surface-soft);
  padding: 12px 14px 4px;
  transition: opacity 0.2s ease, background-color 0.2s ease;
}
.ent-filter-box.ent-filter-disabled {
  opacity: 0.68;
  background: #f3f6f7;
}
.ent-filter-header .form-group {
  margin-bottom: 8px;
}
.ent-filter-header .checkbox {
  margin: 0;
}
.ent-filter-header .checkbox label {
  font-weight: 700;
  color: var(--rve-navy);
  position: relative;
  padding-left: 50px;
  min-height: 24px;
  line-height: 24px;
}
.ent-filter-header #ent_filter_enabled {
  position: absolute;
  opacity: 0;
  pointer-events: none;
}
.ent-filter-header .checkbox label::before {
  content: "";
  position: absolute;
  left: 0;
  top: 1px;
  width: 40px;
  height: 22px;
  border-radius: 999px;
  background: #9aa8ad;
  transition: background-color 0.2s ease;
}
.ent-filter-header .checkbox label::after {
  content: "";
  position: absolute;
  left: 3px;
  top: 4px;
  width: 16px;
  height: 16px;
  border-radius: 50%;
  background: #ffffff;
  box-shadow: 0 1px 3px rgba(16, 32, 51, 0.25);
  transition: transform 0.2s ease;
}
.ent-filter-header .checkbox label:has(#ent_filter_enabled:checked)::before {
  background: var(--rve-teal);
}
.ent-filter-header .checkbox label:has(#ent_filter_enabled:checked)::after {
  transform: translateX(18px);
}
.ent-filter-header .checkbox label:focus-within::before {
  outline: 3px solid rgba(74, 159, 216, 0.28);
  outline-offset: 2px;
}
.ent-filter-fields {
  display: grid;
  grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
  gap: 12px;
}
@media (max-width: 768px) {
  .ent-filter-fields {
    grid-template-columns: 1fr;
    gap: 0;
  }
}
```

- [ ] **Step 4: Run the tests and verify GREEN**

Run:

```bash
Rscript tests/testthat.R
```

Expected: helper and UI layout tests PASS.

- [ ] **Step 5: Parse the modified UI**

Run:

```bash
Rscript -e "parse(file='ui.R'); cat('ui.R parse OK\n')"
```

Expected: `ui.R parse OK`.

- [ ] **Step 6: Commit the layout**

```bash
git add ui.R tests/testthat/test-conservation-filter-ui.R
git commit -m "feat: reorganize Conservation filter controls"
```

### Task 3: Switch State, Modal, and Effective Filtering

**Files:**
- Modify: `tests/testthat/test-conservation-filter-ui.R`
- Modify: `server.R:928-946`
- Modify: `server.R:2598-2635`

- [ ] **Step 1: Add failing server contract tests**

Append:

```r
testthat::test_that("server toggles controls and warns whenever filtering is enabled", {
  server_source <- paste(readLines("server.R", warn = FALSE), collapse = "\n")

  testthat::expect_match(
    server_source,
    'observeEvent\\(input\\$ent_filter_enabled'
  )
  testthat::expect_match(
    server_source,
    'shinyjs::toggleState\\("ent_group_by",\\s*condition = isTRUE\\(input\\$ent_filter_enabled\\)\\)'
  )
  testthat::expect_match(
    server_source,
    'shinyjs::toggleState\\("ent_group",\\s*condition = isTRUE\\(input\\$ent_filter_enabled\\)\\)'
  )
  testthat::expect_match(
    server_source,
    'Real-time calculations using a group filter may take longer'
  )
})

testthat::test_that("Conservation calculations use the effective group", {
  server_source <- paste(readLines("server.R", warn = FALSE), collapse = "\n")

  testthat::expect_match(server_source, 'ent_effective_group <- reactive\\(')
  testthat::expect_match(
    server_source,
    'conservation_effective_group\\(input\\$ent_filter_enabled, input\\$ent_group\\)'
  )
  testthat::expect_match(
    server_source,
    'usage_entropy_data\\([^\\n]*ent_effective_group\\(\\)'
  )
})
```

- [ ] **Step 2: Run the tests and verify RED**

Run:

```bash
Rscript tests/testthat.R
```

Expected: server contract tests FAIL because no state observer, modal, or effective-group reactive exists.

- [ ] **Step 3: Preserve the selected group when choices refresh**

In the `ent_group` choice observer, replace unconditional `selected = "All"` with:

```r
selected_group <- if (
  !is.null(input$ent_group) &&
  input$ent_group %in% c("All", clade_choices)
) {
  input$ent_group
} else {
  "All"
}

updateSelectInput(
  session,
  "ent_group",
  choices = c("All", clade_choices),
  selected = selected_group
)
```

Apply the same preservation logic to both the DuckDB and RDS branches.

- [ ] **Step 4: Implement switch state and modal behavior**

Add:

```r
observe({
  enabled <- isTRUE(input$ent_filter_enabled)
  shinyjs::toggleState("ent_group_by", condition = enabled)
  shinyjs::toggleState("ent_group", condition = enabled)
  shinyjs::toggleClass(
    selector = ".ent-filter-box",
    class = "ent-filter-disabled",
    condition = !enabled
  )
})

observeEvent(input$ent_filter_enabled, {
  if (!isTRUE(input$ent_filter_enabled)) return()

  showModal(modalDialog(
    title = "Group filtering may take longer",
    "Real-time calculations using a group filter may take longer. Please wait while the results are recalculated.",
    easyClose = TRUE,
    footer = modalButton("OK")
  ))
}, ignoreInit = TRUE)
```

- [ ] **Step 5: Route calculations and labels through the effective group**

Add:

```r
ent_effective_group <- reactive({
  conservation_effective_group(input$ent_filter_enabled, input$ent_group)
})
```

Update the title:

```r
output$ent_plot_title <- renderText({
  effective_group <- ent_effective_group()
  clade_text <- if (effective_group == "All") {
    paste("All", input$ent_group_by)
  } else {
    paste(input$ent_group_by, effective_group)
  }
  mode_text <- if (input$variation_type == "AA") "Amino Acid" else "Nucleotide"
  paste(mode_text, "Shannon Entropy Landscape - Subtype", input$global_subtype, "| Gene", input$ent_gene, "|", clade_text)
})
```

In `entropy_site_summary()`, bind `effective_group <- ent_effective_group()`, pass it to `usage_entropy_data()`, and use it instead of `input$ent_group` in the RDS branch filter.

- [ ] **Step 6: Run tests and parse checks**

Run:

```bash
Rscript tests/testthat.R
Rscript -e "parse(file='global.R'); parse(file='ui.R'); parse(file='server.R'); cat('R parse checks OK\n')"
```

Expected: all tests PASS and output ends with `R parse checks OK`.

- [ ] **Step 7: Commit behavior**

```bash
git add server.R tests/testthat/test-conservation-filter-ui.R
git commit -m "feat: enable optional Conservation group filtering"
```

### Task 4: End-to-End Verification

**Files:**
- Modify only if verification exposes a defect: `ui.R`, `server.R`, or tests

- [ ] **Step 1: Start the Shiny app**

Run:

```bash
Rscript -e "shiny::runApp('.', host='127.0.0.1', port=4258, launch.browser=FALSE)"
```

Expected: the app listens on `http://127.0.0.1:4258`.

- [ ] **Step 2: Verify the initial Conservation state in the browser**

Open the app, select `Conservation`, and verify:

- Gene appears to the left of the group-filter box.
- `Filter by group` is off.
- `Group by` and `Group` are disabled and visually muted.
- The plot title identifies all groups.
- The two filter fields sit side by side at desktop width.

- [ ] **Step 3: Verify enabled behavior**

Turn on `Filter by group` and verify:

- The modal appears with the approved performance warning.
- Dismissing with `OK` leaves the filter enabled.
- Both selectors become enabled.
- Selecting a specific group recalculates the entropy plot and updates the title.
- Existing loading feedback remains visible while recalculating.

- [ ] **Step 4: Verify repeated warning and preserved selection**

Turn the switch off, then on again, and verify:

- Disabled state returns to all-group calculations.
- The previously selected `Group by` and `Group` values remain selected.
- The modal appears again on the second enable.

- [ ] **Step 5: Verify responsive stacking**

Set the browser viewport below 768 pixels and verify `Group by` and `Group` stack vertically without overflow.

- [ ] **Step 6: Run final verification**

Run:

```bash
Rscript tests/testthat.R
Rscript -e "parse(file='global.R'); parse(file='ui.R'); parse(file='server.R'); cat('R parse checks OK\n')"
git diff --check
git status --short
```

Expected: tests PASS, parse checks PASS, `git diff --check` reports no errors, and only intentional files are changed.

- [ ] **Step 7: Commit any verification fixes**

If verification required changes:

```bash
git add ui.R server.R R/conservation-filter.R tests/testthat.R tests/testthat
git commit -m "fix: polish Conservation group filter interaction"
```
