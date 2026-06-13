# ==========================================
# 2. UI (User Interface)
# ==========================================

library(shinyjs)
library(bslib)

rve_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#236fa8",
  success = "#168f88",
  info = "#42bfd7",
  danger = "#d9534f",
  bg = "#f4f8fb",
  fg = "#182a3d"
)

rve_card <- function(..., title = NULL, class = "analysis-panel") {
  bslib::card(
    class = paste("rve-bslib-card", class),
    if (!is.null(title)) bslib::card_header(div(class = "analysis-panel-title", title)),
    bslib::card_body(...)
  )
}

rve_control_card <- function(...) {
  bslib::card(
    class = "rve-bslib-card control-panel",
    bslib::card_body(...)
  )
}

rve_stat_card <- function(label, value) {
  bslib::card(
    class = "rve-bslib-card stat-card",
    bslib::card_body(
      h4(label),
      h2(value)
    )
  )
}

ui <- navbarPage(
  id = "main_nav",
  theme = rve_theme,
  
  title = div(
    tags$img(src = "app_icon_round.png", height = "38px", style = "margin-right: 12px; vertical-align: middle;"),
    "RVEAtlas"
  ),
  
  header = tags$head(
    tags$link(rel = "shortcut icon", href = "app_icon_round.png"), 
    use_waiter(), 
    shinyjs::useShinyjs(),
    tags$style(HTML("
      :root {
        --rve-navy: #182a3d;
        --rve-blue: #236fa8;
        --rve-teal: #168f88;
        --rve-ink: #263747;
        --rve-muted: #657482;
        --rve-bg: #f4f8fb;
        --rve-surface: #ffffff;
        --rve-surface-soft: #f8fbfd;
        --rve-border: #dbe5ee;
        --rve-shadow: 0 8px 22px rgba(24, 42, 61, 0.08);
        --rve-shadow-soft: 0 3px 12px rgba(24, 42, 61, 0.06);
        --rve-radius: 8px;
      }
      html, body {
        background: var(--rve-bg);
        color: var(--rve-ink);
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
      }
      body.bslib-page-navbar {
        background:
          linear-gradient(180deg, rgba(232, 243, 251, 0.88) 0, rgba(244, 248, 251, 0) 220px),
          var(--rve-bg);
      }
      @media (min-width: 768px) { 
        .modal-dialog { width: 80vw !important; max-width: 80vw !important; } 
      }
      .jumbotron {
        background: var(--rve-surface);
        padding: 28px 32px;
        border-radius: var(--rve-radius);
        border: 1px solid var(--rve-border);
        box-shadow: var(--rve-shadow);
        margin-top: 20px;
      }
      .jumbotron h1 {
        color: var(--rve-navy) !important;
        font-weight: 800 !important;
        letter-spacing: 0;
        font-size: clamp(38px, 5vw, 64px);
      }
      .jumbotron p {
        color: var(--rve-muted) !important;
      }
      .navbar {
        min-height: 76px;
        background: linear-gradient(135deg, #236fa8 0%, #174f7a 100%);
        border: 0;
        border-bottom: 1px solid var(--rve-border);
        box-shadow: 0 2px 16px rgba(24, 42, 61, 0.06);
        padding: 10px 0;
      }
      .navbar-brand {
        padding-top: 0;
        padding-bottom: 0;
        height: auto;
        color: #ffffff !important;
        font-weight: 800;
        font-size: 26px;
        line-height: 1.1;
        display: flex;
        align-items: center;
        white-space: nowrap;
      }
      .navbar-brand:hover,
      .navbar-brand:focus {
        color: #ffffff !important;
      }
      .navbar-brand img,
      .navbar-brand .fa,
      .navbar-brand svg {
        flex: 0 0 auto;
      }
      .navbar .container-fluid {
        display: flex;
        align-items: center;
        flex-wrap: wrap;
        gap: 14px 18px;
      }
      .navbar-header {
        float: none;
        flex: 0 0 auto;
        order: 1;
        display: flex;
        align-items: center;
        min-height: 56px;
      }
      .navbar-collapse {
        flex: 1 1 500px;
        order: 2;
        padding-right: 0;
        display: flex !important;
        align-items: center;
        justify-content: center;
        min-height: 56px;
      }
      .navbar-nav {
        float: none;
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        justify-content: center;
        gap: 8px;
      }
      .navbar-nav .nav-item {
        margin: 0 2px 6px;
      }
      .navbar .nav-link,
      .navbar .navbar-nav .nav-link {
        color: rgba(255, 255, 255, 0.82) !important;
        font-weight: 700;
        border-radius: 6px;
        font-size: 17px;
      }
      .navbar .nav-link:hover,
      .navbar .nav-link:focus {
        color: #ffffff !important;
        background: rgba(255, 255, 255, 0.12);
      }
      .navbar .nav-link.active,
      .navbar .navbar-nav .active > .nav-link {
        color: #ffffff !important;
        background: rgba(5, 28, 44, 0.28);
        box-shadow: inset 0 -3px 0 #8ee0d4;
      }
      .navbar .container-fluid > .navbar-nav {
        order: 2;
        flex: 1 1 500px;
        min-width: 0;
        justify-content: center;
      }
      .navbar .container-fluid > .navbar-nav > li > a {
        padding: 11px 13px;
        white-space: nowrap;
        color: rgba(255, 255, 255, 0.82) !important;
        font-weight: 700;
        font-size: 17px;
        outline: none;
      }
      .navbar .container-fluid > .navbar-nav > li > a:hover,
      .navbar .container-fluid > .navbar-nav > li > a:focus {
        background: rgba(255, 255, 255, 0.12);
        color: #ffffff !important;
      }
      .navbar .container-fluid > .navbar-nav > .active > a,
      .navbar .container-fluid > .navbar-nav > .active > a:hover,
      .navbar .container-fluid > .navbar-nav > .active > a:focus {
        background: rgba(5, 28, 44, 0.28);
        color: #ffffff !important;
        box-shadow: inset 0 -3px 0 #8ee0d4;
      }
      .navbar .navbar-nav .nav-link.active,
      .navbar .navbar-nav .nav-link.active:hover,
      .navbar .navbar-nav .nav-link.active:focus {
        background: rgba(5, 28, 44, 0.34) !important;
        color: #ffffff !important;
        box-shadow: inset 0 -3px 0 #8ee0d4;
      }
      .tab-content { min-height: calc(100vh - 160px); }
      .tab-content > .tab-pane > .container-fluid {
        padding-top: 22px;
        padding-bottom: 28px;
      }
      .control-panel,
      .well {
        background: var(--rve-surface);
        border: 1px solid var(--rve-border);
        border-radius: var(--rve-radius);
        box-shadow: var(--rve-shadow-soft);
        padding: 18px 20px;
        overflow: visible;
      }
      .well {
        margin-bottom: 18px;
      }
      .year-range-control {
        max-width: 860px;
      }
      .year-range-control .form-group {
        width: 100%;
      }
      .year-range-control .irs {
        width: 100%;
      }
      .analysis-panel {
        background: var(--rve-surface);
        border: 1px solid var(--rve-border);
        border-radius: var(--rve-radius);
        box-shadow: var(--rve-shadow-soft);
        padding: 14px 16px 16px;
        margin: 12px 0 18px;
      }
      .rve-bslib-card {
        border: 1px solid var(--rve-border);
        border-radius: var(--rve-radius);
        box-shadow: var(--rve-shadow-soft);
        background: var(--rve-surface);
        overflow: visible !important;
      }
      .rve-bslib-card.html-fill-container,
      .rve-bslib-card .html-fill-container,
      .rve-bslib-card .html-fill-item {
        overflow: visible !important;
      }
      .rve-bslib-card > .card-body {
        padding: 14px 16px 16px;
        overflow: visible !important;
      }
      .rve-bslib-card > .card-header {
        background: linear-gradient(180deg, #ffffff 0%, #f7fbfe 100%);
        border-bottom: 1px solid var(--rve-border);
        padding: 12px 16px 10px;
        overflow: visible !important;
      }
      .rve-bslib-card .form-group,
      .rve-bslib-card .selectize-control,
      .rve-bslib-card .bootstrap-select {
        overflow: visible !important;
      }
      .selectize-dropdown,
      .bootstrap-select .dropdown-menu {
        z-index: 3000 !important;
      }
      .rve-bslib-card > .card-header .analysis-panel-title {
        margin: 0;
      }
      .analysis-panel-title,
      .dashboard-section-title {
        font-size: 17px;
        font-weight: 750;
        color: var(--rve-navy);
        margin: 2px 0 12px;
      }
      .stat-card {
        min-height: 116px;
        background: linear-gradient(180deg, #ffffff 0%, #f7fbfe 100%);
        border: 1px solid var(--rve-border);
        border-radius: var(--rve-radius);
        box-shadow: var(--rve-shadow-soft);
        padding: 18px 16px;
        text-align: center;
        margin-bottom: 18px;
      }
      .stat-card h4 {
        margin: 0 0 8px;
        color: var(--rve-muted);
        font-size: 14px;
        font-weight: 750;
        text-transform: uppercase;
        letter-spacing: 0.5px;
      }
      .stat-card h2 {
        margin: 0;
        color: var(--rve-blue);
        font-size: 30px;
        font-weight: 800;
      }
      .methods-card {
        max-width: 960px;
        margin: 24px auto;
        padding: 24px 28px;
        background: var(--rve-surface);
        border: 1px solid var(--rve-border);
        border-radius: var(--rve-radius);
        box-shadow: var(--rve-shadow);
      }
      label {
        color: var(--rve-ink);
        font-weight: 700;
      }
      .form-control,
      .selectize-input,
      .bootstrap-select .btn {
        border-color: #cfdce7;
        border-radius: 6px;
      }
      .btn {
        border-radius: 6px;
        font-weight: 700;
      }
      
      /* Prominent navbar pathogen/subtype switch */
      .navbar-pathogen-switch {
        position: relative;
        float: none;
        order: 3;
        margin: 0 12px 0 auto;
        box-sizing: border-box;
        max-width: calc(100vw - 24px);
        background-color: #ffffff;
        background-image: linear-gradient(135deg, #ffffff 0%, #eef7ff 100%);
        padding: 8px 15px;
        border-radius: 34px;
        border: 2px solid var(--rve-blue);
        box-shadow: 0 3px 12px rgba(23, 79, 122, 0.22);
        display: flex;
        align-items: center; 
        gap: 12px;
        min-height: 62px;
        z-index: 10;
      }
      .navbar-pathogen-switch .btn-group-container-sw {
        display: flex;
        align-items: center;
        gap: 10px;
      }
      .switch-label {
        font-weight: 800;
        color: var(--rve-navy);
        font-size: 14px;
        text-transform: uppercase;
        letter-spacing: 1px;
        margin-bottom: 0;
      }
      .navbar-pathogen-switch .bootstrap-select {
        width: auto !important;
        max-width: 100%;
      }
      .navbar-pathogen-switch .bootstrap-select .btn {
        background-color: var(--rve-blue) !important;
        color: white !important;
        border-radius: 22px !important;
        font-weight: bold !important;
        border: none !important;
        min-width: 138px;
        font-size: 18px;
        padding: 9px 17px;
        box-shadow: inset 0 -1px 0 rgba(0,0,0,0.18);
      }
      .navbar-pathogen-switch .bootstrap-select:nth-of-type(2) .btn {
        background-color: var(--rve-teal) !important;
      }
      .navbar-pathogen-switch > .form-group:nth-of-type(2) {
        width: 108px !important;
        flex: 0 1 108px;
      }
      .navbar-pathogen-switch > .form-group:nth-of-type(2) .bootstrap-select .btn {
        min-width: 108px;
      }
      .navbar-pathogen-switch .form-group { margin-bottom: 0 !important; }
      
      /* Custom Radio Button Styling */
      .navbar-pathogen-switch .btn-default {
        background-color: #f1f3f5 !important;
        color: #495057 !important;
        border: 1px solid #ced4da !important;
        font-weight: 700 !important;
        transition: all 0.2s ease-in-out;
        padding: 5px 12px;
      }
      .navbar-pathogen-switch .btn-default.active {
        background-color: #e74c3c !important; /* Bold color for selection */
        color: #ffffff !important;
        border-color: #c0392b !important;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.2) !important;
      }
      .navbar-pathogen-switch .btn-default:hover:not(.active) {
        background-color: #e9ecef !important;
      }
      @media (max-width: 1240px) {
        .navbar .container-fluid { align-items: center; }
        .navbar .container-fluid > .navbar-nav {
          order: 3;
          flex: 1 1 100%;
          justify-content: center;
        }
        .navbar-collapse {
          flex: 1 1 100%;
          order: 3;
          justify-content: center;
          min-height: auto;
        }
        .navbar-pathogen-switch {
          float: none;
          clear: both;
          order: 2;
          margin: 0 12px 0 auto;
          width: auto;
          max-width: 620px;
          justify-content: flex-start;
        }
      }
      @media (max-width: 760px) {
        .navbar .container-fluid {
          justify-content: center;
        }
        .navbar-header,
        .navbar-pathogen-switch {
          margin-left: auto;
          margin-right: auto;
        }
        .navbar-brand {
          font-size: 23px;
        }
      }
      @media (max-width: 640px) {
        .navbar-pathogen-switch {
          display: grid;
          grid-template-columns: auto minmax(0, 1fr);
          width: calc(100vw - 28px);
          max-width: calc(100vw - 28px);
          border-radius: 14px;
          gap: 8px;
          justify-content: stretch;
        }
        .navbar-pathogen-switch .switch-label {
          align-self: center;
          white-space: nowrap;
        }
        .navbar-pathogen-switch .bootstrap-select {
          width: 100% !important;
          min-width: 0 !important;
        }
        .navbar-pathogen-switch > .form-group:nth-of-type(2) {
          width: 100% !important;
          flex: 1 1 auto;
        }
        .navbar-pathogen-switch .bootstrap-select .btn {
          width: 100%;
          min-width: 0;
        }
      }
      .entropy-discovery-panel {
        margin-top: 18px;
        padding: 16px 18px;
        background: var(--rve-surface);
        border: 1px solid var(--rve-border);
        border-radius: var(--rve-radius);
        box-shadow: var(--rve-shadow-soft);
      }
      .entropy-site-section {
        margin-top: 12px;
      }
      .entropy-site-section h4 {
        margin: 0 0 8px 0;
        font-weight: 700;
        color: var(--rve-navy);
      }
      .entropy-site-grid {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
      }
      .ent-site-button {
        border: 1px solid #d4dde6;
        border-radius: 6px;
        background: #f8fafc;
        color: var(--rve-ink);
        padding: 7px 10px;
        min-width: 74px;
        text-align: left;
        line-height: 1.15;
        cursor: pointer;
      }
      .ent-site-button:hover,
      .ent-site-button:focus {
        background: #eef6ff;
        border-color: #3498db;
        outline: none;
      }
      .ent-site-button.high {
        border-color: #e74c3c;
        background: #fff4f2;
      }
      .ent-site-button.mid {
        border-color: #f39c12;
        background: #fff9eb;
      }
      .ent-site-position {
        display: block;
        font-weight: 700;
      }
      .ent-site-entropy {
        display: block;
        font-size: 0.86em;
        color: #5f6c7b;
      }
      .entropy-empty {
        color: var(--rve-muted);
        font-style: italic;
      }
      .ent-plot-loading-wrap {
        position: relative;
      }
      .ent-plot-loading-wrap .shiny-bound-output.recalculating {
        opacity: 0.35;
      }
      .ent-plot-loading-wrap .ent_plot_recalculating {
        display: none;
        position: absolute;
        inset: 0;
        z-index: 10;
        align-items: center;
        justify-content: center;
        flex-direction: column;
        gap: 10px;
        background: rgba(255, 255, 255, 0.72);
        color: var(--rve-navy);
        font-weight: 700;
        pointer-events: none;
      }
      .ent-plot-loading-wrap:has(.shiny-bound-output.recalculating) .ent_plot_recalculating,
      .ent-plot-loading-wrap.ent-loading-active .ent_plot_recalculating {
        display: flex;
      }
      .ent-loading-spinner {
        width: 34px;
        height: 34px;
        border-radius: 50%;
        border: 4px solid #d7e6f5;
        border-top-color: #3498db;
        animation: entSpin 0.8s linear infinite;
      }
      @keyframes entSpin {
        to { transform: rotate(360deg); }
      }
      .dataset-breakdown-plot-wrap,
      .dataset-breakdown-plot-wrap .shiny-plot-output,
      .dataset-breakdown-plot-wrap .plotly,
      .dataset-breakdown-plot-wrap .plot-container,
      .dataset-breakdown-plot-wrap .svg-container {
        width: 100% !important;
        max-width: 100% !important;
      }
      .info-markdown table {
        width: 100%;
        table-layout: fixed;
        margin: 1.25rem 0;
        border-collapse: separate;
        border-spacing: 0;
        border: 1px solid #dfe6e9;
        border-radius: 10px;
        overflow: hidden;
      }
      .info-markdown thead th {
        background-color: #f4f7f9;
        color: #2c3e50;
        font-weight: 700;
        text-align: left;
        padding: 14px 20px;
        border-bottom: 2px solid #dfe6e9;
        white-space: nowrap;
      }
      .info-markdown tbody td {
        padding: 12px 20px;
        border-bottom: 1px solid #ecf0f1;
        vertical-align: top;
      }
      .info-markdown tbody tr:last-child td {
        border-bottom: none;
      }
      .info-markdown th + th,
      .info-markdown td + td {
        border-left: 1px solid #ecf0f1;
      }
      .info-markdown code {
        white-space: normal;
        overflow-wrap: anywhere;
        word-break: break-word;
      }
      .summary-card {
        background: var(--rve-surface);
        border: 1px solid var(--rve-border);
        border-radius: var(--rve-radius);
        min-height: 120px;
        padding: 18px 14px;
        text-align: center;
        box-shadow: 0 1px 3px rgba(32,56,85,0.08);
      }
      .summary-card-title {
        font-size: 17px;
        font-weight: 700;
        color: var(--rve-navy);
        margin-bottom: 10px;
      }
      .summary-card-value {
        font-size: 30px;
        font-weight: 800;
        color: var(--rve-blue);
        line-height: 1.15;
      }
      .gc-summary-card {
        height: 150px;
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
      }
      .gc-summary-card-note {
        min-height: 20px;
        margin: 8px 0 0 0;
        color: #5f6c7b;
        font-size: 13px;
      }
    ")),
    tags$script(HTML("
      window.resizeDatasetBreakdownPlot = function() {
        var el = document.getElementById('stats_clade_plot');
        if (!el || !window.Plotly) return;
        setTimeout(function() {
          Plotly.Plots.resize(el);
          window.dispatchEvent(new Event('resize'));
        }, 80);
        setTimeout(function() {
          Plotly.Plots.resize(el);
        }, 350);
      };
      Shiny.addCustomMessageHandler('resize_dataset_breakdown_plot', function(message) {
        window.resizeDatasetBreakdownPlot();
      });
      window.mountNavbarPathogenSwitch = function() {
        var switchEl = document.querySelector('.navbar-pathogen-switch');
        var navInner = document.querySelector('.navbar .container-fluid');
        if (switchEl && navInner && !navInner.contains(switchEl)) {
          var collapse = navInner.querySelector('.navbar-collapse');
          if (collapse) {
            navInner.insertBefore(switchEl, collapse);
          } else {
            navInner.appendChild(switchEl);
          }
        }
      };
      document.addEventListener('DOMContentLoaded', window.mountNavbarPathogenSwitch);
      $(document).on('shiny:connected', window.mountNavbarPathogenSwitch);
      setTimeout(window.mountNavbarPathogenSwitch, 100);
      setTimeout(window.mountNavbarPathogenSwitch, 600);
      window.showEntropyPlotLoading = function() {
        var wrap = document.querySelector('.ent-plot-loading-wrap');
        if (!wrap) return;
        wrap.classList.add('ent-loading-active');
      };
      window.hideEntropyPlotLoading = function() {
        var wrap = document.querySelector('.ent-plot-loading-wrap');
        if (!wrap) return;
        wrap.classList.remove('ent-loading-active');
      };
      $(document).on('plotly_afterplot', '#ent_plot', function() {
        window.hideEntropyPlotLoading();
      });
      $(document).on('shiny:value shiny:outputinvalidated', function(event) {
        if (event.target && event.target.id === 'ent_plot') {
          window.showEntropyPlotLoading();
        }
      });
      $(document).on('click', '.navbar-nav a', function() {
        if ($(this).text().trim() === 'Conservation' || $(this).text().trim() === 'Conservation (Entropy)') {
          window.showEntropyPlotLoading();
        }
      });
      $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"]', function(e) {
        var tabLabel = $(e.target).text().trim();
        if (tabLabel === 'Dataset Insights') {
          window.resizeDatasetBreakdownPlot();
        }
        if (tabLabel === 'Conservation' || tabLabel === 'Conservation (Entropy)') {
          window.showEntropyPlotLoading();
          setTimeout(function() {
            var plot = document.getElementById('ent_plot');
            if (plot && plot.querySelector('.plotly')) {
              window.hideEntropyPlotLoading();
            }
          }, 900);
        }
      });
    ")),
    div(class = "navbar-form navbar-right navbar-pathogen-switch",
        span(class = "switch-label", "Pathogen"),
        pickerInput(
          inputId = "active_pathogen",
          label = NULL,
          choices = pathogen_choices(),
          selected = "FLU",
          width = "140px"
        ),
        span(class = "switch-label", "Subtype"),
        pickerInput(
          inputId = "global_subtype",
          label = NULL,
          choices = metadata_groups,
          selected = metadata_groups[1],
          width = "108px"
        ),
        div(style = "width: 1px; height: 25px; background: #ced4da;display:None"),
        div(class = "btn-group-container-sw", style = "display:None",
            span(class = "switch-label", "Mode"),
            radioGroupButtons(
              inputId = "variation_type",
              label = NULL,
              choices = c("AA", "NT"),
              selected = "AA",
              status = "default",
              size = "sm"
            )
        )
    )
  ),
  
  footer = tags$footer(
    style = "text-align: center; padding: 15px; background-color: #f8f9fa; border-top: 1px solid #e7e7e7; color: #6c757d; margin-top: 30px; width: 100%;",
    HTML(paste0("&copy; ", format(Sys.Date(), "%Y"), "Center for Applied Bioinformatics | St. Jude Research. All rights reserved."))
  ),
  
  # ---------------------------------------------------------
  # MEMORY MONITOR & CONTROL WIDGET
  # ---------------------------------------------------------
  tags$div(
    style = "position: fixed; bottom: 15px; left: 15px; z-index: 9999; background-color: rgba(255,255,255,0.95); padding: 6px 15px; border-radius: 30px; box-shadow: 0 4px 10px rgba(0,0,0,0.15); font-size: 13px; display: flex; align-items: center; gap: 12px; border: 1px solid #dee2e6;",
    tags$strong(icon("memory"), "RAM:"),
    textOutput("mem_usage", inline = TRUE),
    actionButton("free_mem", "Clear Cache", icon = icon("broom"), class = "btn-xs btn-danger", style = "border-radius: 20px; padding: 2px 10px; font-size: 12px; font-weight: bold;")
  ),
  
  # ---------------------------------------------------------
  # TAB 0: HOME
  # ---------------------------------------------------------
  tabPanel("Home", value = "home",
           fluidPage(
             div(class = "jumbotron",
                 h1("Welcome to the RVEAtlas Explorer", style = "color: #2c3e50; font-weight: bold;"),
                 p("A high-resolution visualization tool for analyzing Respiratory Virus genetic diversity across multiple subtypes and lineages.", style = "font-size: 1.2em; color: #7f8c8d;"),
                 hr(),
                 div(style = "text-align: center; margin-top: 20px; margin-bottom: 30px;",
                     uiOutput("home_welcome_banner")
                 ),
                 h3("How to Use This App:", style = "color: #2980b9;"),
                 fluidRow(
                   column(4, h4(icon("chart-bar"), " Single Position Explorer"), p("Dive deep into the amino acid or nucleotide distribution of any specific position within an viral gene.")),
                   column(4, h4(icon("not-equal"), " Pairwise Comparison"), p("Instantly identify robust, fixed differences between any two genetic clades across all genes.")),
                   column(4, h4(icon("globe"), " Gene-Wide Landscapes"), p("Explore whole-gene visualizations including Entropy conservation plots and genetic clade dynamics"))
                 )
             )
           )
  ),
  
  # TAB 1: DATASET STATS (World Map + Static Plots)
  tabPanel("Dataset Insights", value = "dataset_insights",
           fluidPage(
             fluidRow(
               column(4, rve_stat_card("Total Sequences", textOutput("total_seqs"))),
               column(4, rve_stat_card("Countries Represented", textOutput("total_countries"))),
               column(4, rve_stat_card("Time Span", textOutput("time_range")))
             ),
             
             rve_control_card(
               fluidRow(
                 # column(2, selectInput("map_geo_level", "Grouping:", choices = c("Region", "Country"))),
                 # column(2, selectInput("map_clade_type", "Pie Data:", choices = c("HA-Clade" = "clade", "NA-Clade" = "G_clade"))),
                 # column(3, selectInput("map_year", "Select Year:", choices = c("All", metadata_years))),
                 column(12, class = "year-range-control",
                        sliderInput("stats_year_range", "Plot Year Range:",
                                    min = 1918, 
                                    max = as.numeric(format(Sys.Date(), "%Y")), 
                                    value = c(1968, as.numeric(format(Sys.Date(), "%Y"))), 
                                    sep = "", step = 1))
               ),
               helpText("Note: 'Plot Year Range' affects the charts below.")
             ),

             # fluidRow(
             #   column(12, 
             #          h4("Global Clade Distribution", style="font-weight: bold; color: #2c3e50;"),
             #          leafletOutput("world_map", height = "500px"),
             #          hr())
             # ),

             
             fluidRow(
               column(6, 
                      rve_card(
                          title = "Sequencing Over Time (Seasonality)",
                          withWaiter(plotlyOutput("stats_time_plot", height = "400px"))
                      )
               ),
               column(6, 
                      rve_card(
                          title = "Regional Breakdown",
                          withWaiter(plotlyOutput("stats_geo_plot", height = "400px"))
                      )
               )
             ),
             
             fluidRow(
               column(12, hr()),
               column(3, selectInput("clade_plot_subtype", "Filter Subtype:", 
                                     choices = metadata_groups, 
                                     selected = metadata_groups[1])),
               column(3, selectInput("clade_plot_fill", "Sub-Category (Color):", choices = NULL)),
               column(3, selectInput("clade_plot_palette", "Color Palette:", 
                                     choices = c("Viridis" = "viridis", "Plasma" = "plasma", 
                                                 "Magma" = "magma", "Inferno" = "inferno", 
                                                 "Cividis" = "cividis", "Turbo" = "turbo", "Rainbow" = "rainbow"), selected = "turbo")),
               column(3,
                      selectInput("clade_plot_time_scale", "X Axis:",
                                  choices = c("Year" = "Year", "Year-Month" = "YearMonth"),
                                  selected = "Year",
                                  selectize = FALSE)),
               column(12,
                      rve_card(
                          title = textOutput("stats_clade_plot_title", inline = TRUE),
                          class = "analysis-panel dataset-breakdown-plot-wrap",
                          withWaiter(plotlyOutput("stats_clade_plot", height = "500px", width = "100%"))
                      )
               )
             )
           )
  ),
  
  # ---------------------------------------------------------
  # MACRO-LEVEL DROPDOWN MENU
  # ---------------------------------------------------------
  # navbarMenu("Gene-Wide Landscapes",
             
             tabPanel("Conservation", value = "conservation",
                      fluidPage(
                        rve_control_card(
                          fluidRow(
                            column(3, helpText("Calculates Shannon Entropy to identify highly conserved valleys and hypervariable peaks across the entire gene. Subtype is controlled globally.")),
                            column(3, selectInput("ent_group_by", "Group by:", choices = NULL)),
                            column(3, selectInput("ent_gene", "Gene:", choices = NULL)),
                            column(3, selectInput("ent_group", "Group:", choices = NULL))
                          ),
                          fluidRow(
                            column(3, sliderInput("ent_min_seqs", "Min Sequences:", min = 0, max = 1000, value = 10, step = 10)),
                            column(3, sliderInput("ent_font_size", "Plot Font Size:", min = 10, max = 24, value = 14, step = 1)),
                            column(3, radioButtons("ent_plot_format", "Format:", choices = c("PNG", "PDF"), inline = TRUE)),
                            column(3, downloadButton("downloadEntPlot", "Download Plot", class = "btn-info", style="margin-top: 25px; width: 100%;"))
                          )
                        ),
                        rve_card(
                            title = textOutput("ent_plot_title"),
                            div(class = "ent-plot-loading-wrap",
                                withWaiter(plotlyOutput("ent_plot", height = "450px")),
                                div(class = "ent_plot_recalculating",
                                    div(class = "ent-loading-spinner"),
                                    div("Updating entropy plot...")
                                )
                            )
                        ),
                        div(class = "entropy-discovery-panel",
                            h3("Variable Sites", style = "margin-top: 0; font-weight: bold; color: #2c3e50;"),
                            p("Click a position to inspect its full distribution in the Single Position Explorer.", style = "color: #5f6c7b; margin-bottom: 8px;"),
                            uiOutput("ent_variable_sites")
                        )
                      )
             ),
             # 
             # tabPanel("Mutation Tracker (Lollipop)",
             #          fluidPage(
             #            wellPanel(
             #              fluidRow(
             #                column(3, helpText("Visualize fixed amino acid mutations in a Target Group compared to a Reference Group. Subtype is controlled globally.")),
             #                column(3, selectInput("lol_group_by", "Group by:", choices = NULL)),
             #                column(3, selectInput("lol_gene", "Gene:", choices = NULL)),
             #                column(3, numericInput("lol_min_freq", "Min Dominant Freq (%):", value = 90.0, min = 50.0, max = 100.0))
             #              ),
             #              fluidRow(
             #                column(3, selectInput("lol_ref_group", "Reference Group:", choices = NULL)),
             #                column(3, selectInput("lol_tar_group", "Target Group:", choices = NULL)),
             #                column(2, sliderInput("lol_font_size", "Font Size:", min = 10, max = 24, value = 14, step = 1)),
             #                column(2, radioButtons("lol_plot_format", "Format:", choices = c("PNG", "PDF"), inline = TRUE)),
             #                column(2, downloadButton("downloadLolPlot", "Download Plot", class = "btn-info", style="margin-top: 25px; width: 100%;"))
             #              )
             #            ),
             #            h3(textOutput("lol_plot_title")),
             #            withWaiter(plotlyOutput("lol_plot", height = "550px")) 
             #          )
             # ),
             
            #  # TAB 5: CONSENSUS MSA (FULL-WIDTH LAYOUT)
            #  tabPanel("Consensus MSA Map",
            #           fluidPage(
            #             wellPanel(
            #               fluidRow(
            #                 column(3, helpText("Interactive Multiple Sequence Alignment. Subtype is controlled globally.")),
            #                 column(3, selectInput("heat_group_by", "Group by:", choices = NULL)),
            #                 column(3, selectInput("heat_gene", "Gene:", choices = NULL)),
            #                 column(3, div(style = "margin-top: 25px;", checkboxInput("show_mut_only", "Show Mutations Only", value = FALSE)))
            #               )
            #             ),
            #             fluidRow(
            #               column(12,
            #                      h3(textOutput("heat_plot_title")),
            #                      withWaiter(uiOutput("msa_dynamic_container")) 
            #               )
            #             )
            #           )
            #  )
  # ),

  # ---------------------------------------------------------
  # TAB 1B: GENETIC CLADE
  # ---------------------------------------------------------
  tabPanel("Genetic Clade", value = "genetic_clade",
           fluidPage(
             rve_control_card(
               fluidRow(
                 column(4, selectInput("gc_annotation", "Clade annotation:", choices = NULL)),
                 column(
                   8,
                   selectizeInput(
                     "gc_clade",
                     "Clade / group:",
                     choices = NULL,
                     selected = NULL,
                     options = list(
                       placeholder = "Type to search or choose a clade",
                       maxOptions = 500,
                       maxItems = 1,
                       create = FALSE,
                       searchField = c("label", "value"),
                       openOnFocus = TRUE,
                       closeAfterSelect = TRUE
                     )
                   )
                 )
               )
             ),
             uiOutput("gc_status_notice"),
             uiOutput("gc_summary_cards"),
             fluidRow(
               column(
                 12,
                 rve_card(
                     title = "Monthly Clade Prevalence",
                     withWaiter(plotlyOutput("gc_prevalence_plot", height = "430px"))
                 )
               )
             ),
             fluidRow(
               column(4, rve_card(title = "Top Countries", withWaiter(plotlyOutput("gc_country_plot", height = "300px")))),
               column(4, rve_card(title = "Top Regions", withWaiter(plotlyOutput("gc_region_plot", height = "300px")))),
               column(4, rve_card(title = "Top Hosts", withWaiter(plotlyOutput("gc_host_plot", height = "300px"))))
             ),
             fluidRow(
               column(
                 12,
                 rve_card(
                     title = "Monthly Detail",
                     DTOutput("gc_monthly_table")
                 )
               )
             )
           )
  ),
  
  # ---------------------------------------------------------
  # TAB 1: SINGLE POSITION
  # ---------------------------------------------------------
  tabPanel("Single Position", value = "single_position",
           fluidPage(
             rve_control_card(
               fluidRow(
                 column(3,
                        h5("Setting", style="font-weight: bold; color: #2980b9;"),
                        selectInput("sp_group_by", "Group by:", choices = NULL),
                        uiOutput("sp_year_month_range_ui"),
                        checkboxInput("sp_hide_empty_years", "Hide years without records (when Group by Year)", value = TRUE),
                        checkboxInput("sp_show_counts", "Show raw counts instead of percentage", value = FALSE)
                 ),
                 column(2,
                        div(id = "sp_quick_access_section",
                            h5("Quick Access", style="font-weight: bold; color: #2980b9;"),
                            selectInput("sp_quick_visit", "Jump to Key Position:", 
                                        choices = if(nrow(important_pos_df) > 0) c("Manual Selection" = "None", setNames(1:nrow(important_pos_df), important_pos_df$label)) else c("Manual Selection" = "None"))
                        ),
                        sliderInput("sp_min_seqs", "Min Seqs:", min = 1, max = 500, value = 5)
                 ),
                 column(5,
                        h5("Precise Access", style="font-weight: bold; color: #2980b9;"),
                        fluidRow(
                          column(12, selectInput("sp_gene", "Gene:", choices = NULL))
                        ),
                        fluidRow(
                          column(12,
                                 uiOutput("sp_range_label"),
                                 div(class = "sp-position-row",
                                     tags$style(HTML("
                                       .sp-position-row {
                                         width: 100%;
                                         max-width: 560px;
                                         margin-bottom: 15px;
                                       }
                                       #sp_position {
                                         width: 100% !important;
                                         margin-bottom: 0px !important;
                                       }
                                       #sp_position .selectize-control,
                                       #sp_position .selectize-input {
                                         width: 100% !important;
                                       }
                                       #sp_position input {
                                         height: 34px;
                                         padding-top: 6px;
                                         padding-bottom: 6px;
                                       }
                                       .sp-numbering-hint {
                                         min-height: 22px;
                                         margin-top: 6px;
                                       }
                                     ")),
                                     selectizeInput("sp_position", label = NULL, choices = NULL, selected = NULL,
                                                    width = "100%",
                                                    options = list(
                                                      placeholder = "Type/select position",
                                                      create = FALSE,
                                                      maxOptions = 1000
                                                    )),
                                     div(class = "sp-numbering-hint", uiOutput("sp_numbering_label"))
                                 )
                          )
                        )
                 ),
                 column(2,
                        h5("Export", style="font-weight: bold; color: #2980b9;"),
                        sliderInput("sp_font_size", "Plot Font Size:", min = 10, max = 24, value = 14),
                        radioButtons("sp_plot_format", "Download Format:", choices = c("PDF", "PNG"), inline = TRUE),
                        downloadButton("downloadSpPlot", "Download Plot", class = "btn-info", style="width: 100%; margin-bottom: 6px;"),
                        downloadButton("downloadSpPositionExcel", "Download Excel Matrix", class = "btn-success", style="width: 100%;")
                 )
               )
             ),
             fluidRow(
               column(12,
                      rve_card(
                          uiOutput("sp_position_details"),
                          uiOutput("sp_position_count_info"),
                          div(style = "margin-bottom: 28px;",
                              plotlyOutput("sp_overall_aa_bar", height = "120px")
                          ),
                          tags$hr(style = "border-top: 1px solid #d8e3ec; margin: 4px 0 28px 0;"),
                          withWaiter(plotlyOutput("sp_aa_plot", height = "500px")),
                          DTOutput("sp_aa_table")
                      )
               )
             )
           )
  ),
  
  # ---------------------------------------------------------
  # TAB 2: PAIRWISE COMPARISON
  # ---------------------------------------------------------
  tabPanel("Pairwise Comparison", value = "pairwise_comparison",
           fluidPage(
             rve_control_card(
               fluidRow(
                 column(3,
                        selectInput("pw_group_by", "Group by:", choices = NULL),
                        checkboxInput("pw_hide_empty_years", "Hide years without records", value = TRUE)
                 ),
                 column(3,
                        selectInput("pw_clade1", "Group 1:", choices = NULL),
                        selectInput("pw_clade2", "Group 2:", choices = NULL)
                 ),
                 column(3,
                        numericInput("pw_min_freq", "Min Dominant Freq (%):", value = 90.0, min = 50.0, max = 100.0, step = 1.0)
                 ),
                 column(3,
                        br(),
                        downloadButton("downloadPairwiseCSV", "Download Table (CSV)", class = "btn-primary", style="margin-bottom: 5px; width: 100%;"),
                        downloadButton("downloadPairwiseExcel", "Download Excel Matrices", class = "btn-success", style="width: 100%;")
                 )
               )
             ),
             fluidRow(
               column(12,
                      rve_card(
                          title = "Cross-Gene Pairwise Differences",
                          p("Click on any highlighted Position to view the full amino acid distribution for that specific site."),
                          DTOutput("pw_diff_table")
                      )
               )
             )
           )
  ),
  # ---------------------------------------------------------
  # TAB 3: METHODS / INFO / UPDATE LOG
  # ---------------------------------------------------------
  tabPanel("Methods & Info", value = "methods_info",
           fluidPage(
             div(
               class = "info-markdown methods-card",
               uiOutput("app_info_markdown")
             )
           )
  )
)
