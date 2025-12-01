# =========================================================
# Shiny App: Auto-clean → Auto-impute → Summaries & Plots
# =========================================================

# Increase upload limit (e.g., 200 MB)
options(shiny.maxRequestSize = 200 * 1024^2)

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)     # fast read & write
  library(tools)     # file_path_sans_ext
  library(GGally)
  library(scales)
  library(tidyr)
  library(e1071)
  library(purrr)
  library(forcats)
  library(stringr)
})

# --- Load your pipeline functions ---
# Put your real file names here or comment these out if functions are already in memory
source("imputation.R")   # should define: auto_clean_mixed_data(), auto_impute_select()
source("viz.R")          # should define: summarize_distributions_comparison(), plot_distributions_comparison(), correlation_compare_plots()

# ---- UI ----
ui <- fluidPage(
  titlePanel("Auto Imputation Pipeline — Upload → Clean → Impute → Compare"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "file", "Upload CSV",
        accept = c(".csv", "text/csv", "text/comma-separated-values")
      ),
      checkboxInput("header", "Header", TRUE),
      selectInput("sep", "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ","),
      selectInput("quote", "Quote", choices = c(None = "", "Double quote" = "\"", "Single quote" = "'"), selected = "\""),
      tags$hr(),
      h5("Imputation Settings"),
      numericInput("m", "m (number of imputations)", value = 5, min = 1, step = 1),
      numericInput("max_iter", "Max iterations", value = 10, min = 1, step = 1),
      numericInput("eval_frac", "Evaluation fraction", value = 0.10, min = 0.05, max = 0.9, step = 0.01),
      numericInput("eval_reps", "Evaluation repetitions", value = 3, min = 1, step = 1),
      textInput("methods", "Methods (comma-separated)", value = "pmm,linreg,linreg_boot,bayes,cart,rf"),
      tags$hr(),
      actionButton("run", "Run Pipeline", class = "btn btn-primary"),
      tags$hr(),
      downloadButton("dl_best_imputed_csv", "Download Imputed CSV"),
      br(), br(),
      downloadButton("dl_missing_comp_csv", "Download: Missing-Features Comparison"),
      br(), br(),
      downloadButton("dl_nonmissing_csv", "Download: Non-Missing Summary")
    ),
    
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Preview",
                           h4("Raw Uploaded Preview"),
                           DTOutput("preview_head")
                  ),
                  
                  tabPanel("Cleaning",
                           h4("auto_clean_mixed_data() — Output Snapshot"),
                           verbatimTextOutput("clean_dim"),
                           DTOutput("clean_head")
                  ),
                  
                  tabPanel("Imputation Results",
                           fluidRow(
                             column(6,
                                    h4("Best Method"),
                                    verbatimTextOutput("best_method"),
                                    h5("Leaderboard"),
                                    DTOutput("leaderboard")
                             ),
                             column(6,
                                    h4("Best Imputed — Head"),
                                    DTOutput("best_imputed_head")
                             )
                           )
                  ),
                  
                  tabPanel("Distributions & Summaries",
                           h4("Comparison for Missing Features (Before vs After)"),
                           DTOutput("tbl_missing_comparison"),
                           tags$hr(),
                           h4("Summary for Non-Missing Features"),
                           DTOutput("tbl_nonmissing_summary")
                  ),
                  tabsetPanel(
                    id = "plot_tabs",
                    tabPanel("Missing-feature comparisons", uiOutput("tabs_missing")),
                    tabPanel("Non-missing summaries", uiOutput("tabs_nonmissing"))
                  )
                  ,
                  
                  tabPanel("Session Log",
                           verbatimTextOutput("log")
                  )
                  
      )
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  # Reactive container for all state
  rv <- reactiveValues(
    raw = NULL,
    df = NULL,
    out = NULL,
    best_method = NULL,
    best_imputed_df = NULL,
    result = NULL,
    plots = NULL,
    cc = NULL,
    log = character()
  )
  
  # Helper to append to Session Log
  append_log <- function(...) {
    rv$log <- c(rv$log, paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = " ")))
  }
  
  # ---- Preview uploaded CSV ----
  observeEvent(input$file, {
    req(input$file)
    append_log("Reading CSV:", input$file$name)
    # Use readr::read_csv for speed and auto types; fall back to base if needed
    tryCatch({
      delim <- input$sep
      quo <- if (nzchar(input$quote)) input$quote else "\""
      # readr read_delim handles header automatically
      raw <- readr::read_delim(
        file = input$file$datapath,
        delim = delim,
        quote = quo,
        col_names = input$header,
        show_col_types = FALSE,
        progress = FALSE
      )
      rv$raw <- as.data.frame(raw)
      append_log("Loaded rows:", nrow(rv$raw), "cols:", ncol(rv$raw))
    }, error = function(e) {
      append_log("ERROR reading CSV:", e$message)
      showNotification(paste("Error reading CSV:", e$message), type = "error", duration = 8)
    })
  })
  
  output$preview_head <- renderDT({
    req(rv$raw)
    datatable(head(rv$raw, 20), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ---- Run Entire Pipeline ----
  observeEvent(input$run, {
    req(rv$raw)
    set.seed(123)
    
    # Parse methods safely
    methods <- strsplit(input$methods, "\\s*,\\s*")[[1]]
    methods <- methods[nzchar(methods)]
    if (length(methods) == 0) {
      showNotification("No methods specified. Using defaults: pmm, linreg, linreg_boot, bayes, cart, rf", type = "warning")
      methods <- c("pmm", "linreg", "linreg_boot", "bayes", "cart", "rf")
    }
    
    withProgress(message = "Running pipeline…", value = 0, {
      incProgress(0.05, "Cleaning data")
      tryCatch({
        rv$df <- auto_clean_mixed_data(rv$raw)
        append_log("auto_clean_mixed_data() OK: rows", nrow(rv$df), "cols", ncol(rv$df))
      }, error = function(e) {
        append_log("ERROR auto_clean_mixed_data():", e$message)
        showNotification(paste("auto_clean_mixed_data() failed:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      
      incProgress(0.35, "Imputing (auto_impute_select)")
      tryCatch({
        rv$out <- auto_impute_select(
          data      = rv$df,
          methods   = methods,
          m         = input$m,
          max_iter  = input$max_iter,
          eval_frac = input$eval_frac,
          eval_reps = input$eval_reps,
          verbose   = TRUE
        )
        append_log("auto_impute_select() finished")
      }, error = function(e) {
        append_log("ERROR auto_impute_select():", e$message)
        showNotification(paste("auto_impute_select() failed:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      
      # Extract best method + best imputed df robustly
      incProgress(0.55, "Selecting best imputed dataset")
      rv$best_method <- tryCatch(rv$out$best_method, error = function(e) NA_character_)
      best_imputed_raw <- tryCatch(rv$out$best_imputed, error = function(e) NULL)
      
      # Cases:
      # (a) best_imputed is already a data.frame
      # (b) best_imputed is a list of m imputations -> choose one deterministically (e.g., first)
      # (c) best_imputed is a named list by methods -> choose by best_method then first element
      rv$best_imputed_df <- NULL
      try({
        if (is.data.frame(best_imputed_raw)) {
          rv$best_imputed_df <- best_imputed_raw
        } else if (is.list(best_imputed_raw)) {
          # If named by method
          if (!is.null(names(best_imputed_raw)) && !is.na(rv$best_method) && rv$best_method %in% names(best_imputed_raw)) {
            x <- best_imputed_raw[[rv$best_method]]
            rv$best_imputed_df <- if (is.list(x) && !is.data.frame(x)) x[[1]] else x
          } else {
            # Otherwise just take the first data frame we find
            first_df <- NULL
            for (x in best_imputed_raw) {
              if (is.data.frame(x)) { first_df <- x; break }
              if (is.list(x) && length(x) && is.data.frame(x[[1]])) { first_df <- x[[1]]; break }
            }
            rv$best_imputed_df <- first_df
          }
        }
      })
      if (is.null(rv$best_imputed_df)) {
        append_log("WARNING: Could not resolve a best imputed data.frame; falling back to raw df")
        rv$best_imputed_df <- rv$df
      } else {
        append_log("Selected best_imputed_df: rows", nrow(rv$best_imputed_df), "cols", ncol(rv$best_imputed_df))
      }
      
      # Summaries
      incProgress(0.75, "Summarizing distributions")
      tryCatch({
        rv$result <- summarize_distributions_comparison(
          before_data = rv$df,
          after_data  = rv$best_imputed_df
        )
        append_log("summarize_distributions_comparison() OK")
      }, error = function(e) {
        append_log("ERROR summarize_distributions_comparison():", e$message)
        showNotification(paste("summarize_distributions_comparison() failed:", e$message), type = "error", duration = 10)
        rv$result <- NULL
      })
      
      # Plots
      incProgress(0.9, "Generating plots")
      tryCatch({
        rv$plots <- plot_distributions_comparison(rv$df, rv$best_imputed_df, cat_top_k = 8)
        append_log("plot_distributions_comparison() OK")
      }, error = function(e) {
        append_log("ERROR plot_distributions_comparison():", e$message)
        showNotification(paste("plot_distributions_comparison() failed:", e$message), type = "error", duration = 10)
        rv$plots <- NULL
      })
      
      tryCatch({
        rv$cc <- correlation_compare_plots(rv$df, rv$best_imputed_df, method = "pearson", reorder = TRUE)
        append_log("correlation_compare_plots() OK")
      }, error = function(e) {
        append_log("ERROR correlation_compare_plots():", e$message)
        showNotification(paste("correlation_compare_plots() failed:", e$message), type = "error", duration = 10)
        rv$cc <- NULL
      })
      
      incProgress(1, "Done")
      updateTabsetPanel(session, "tabs", selected = "Imputation Results")
    })
  })
  
  # ---- Cleaning tab ----
  output$clean_dim <- renderText({
    req(rv$df)
    paste0("Dimensions: ", nrow(rv$df), " × ", ncol(rv$df))
  })
  output$clean_head <- renderDT({
    req(rv$df)
    datatable(head(rv$df, 20), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ---- Imputation Results tab ----
  output$best_method <- renderText({
    req(rv$out)
    paste("Best method:", rv$best_method %||% "Unknown")
  })
  output$leaderboard <- renderDT({
    req(rv$out)
    lb <- tryCatch(rv$out$leaderboard, error = function(e) NULL)
    req(lb)
    datatable(lb, options = list(scrollX = TRUE, pageLength = 10))
  })
  output$best_imputed_head <- renderDT({
    req(rv$best_imputed_df)
    datatable(head(rv$best_imputed_df, 20), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ---- Distributions & Summaries ----
  output$tbl_missing_comparison <- renderDT({
    req(rv$result)
    x <- tryCatch(rv$result$comparison_missing_features, error = function(e) NULL)
    req(x)
    datatable(as.data.frame(x), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$tbl_nonmissing_summary <- renderDT({
    req(rv$result)
    x <- tryCatch(rv$result$summary_nonmissing_features, error = function(e) NULL)
    req(x)
    datatable(as.data.frame(x), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ---- Plots ----
  
  output$p_distributions <- renderPlot({
    req(rv$plots)
    # rv$plots can be a patchwork object or list of ggplots; print whatever it is
    if (inherits(rv$plots, "ggplot") || inherits(rv$plots, "patchwork")) {
      print(rv$plots)
    } else if (is.list(rv$plots)) {
      # if list of ggplots, combine first few
      plist <- rv$plots[1:min(6, length(rv$plots))]
      # Try to wrap using patchwork if possible
      tryCatch({
        wrap <- Reduce(`+`, lapply(plist, function(p) p + theme(plot.margin = margin(5,5,5,5))))
        print(wrap)
      }, error = function(e) {
        # fallback: just print the first
        print(plist[[1]])
      })
    }
  })
  
  output$p_corr <- renderPlot({
    req(rv$cc)
    if (!is.null(rv$cc$plots$combined)) {
      print(rv$cc$plots$combined)
    } else if (!is.null(rv$cc$combined)) {
      print(rv$cc$combined)
    }
  })
  
  # ---- Downloads ----
  output$dl_best_imputed_csv <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$file$name)) file_path_sans_ext(input$file$name) else "imputed"
      paste0(base, "_best_imputed.csv")
    },
    content = function(file) {
      req(rv$best_imputed_df)
      readr::write_csv(rv$best_imputed_df, file)
    }
  )
  
  output$dl_missing_comp_csv <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$file$name)) file_path_sans_ext(input$file$name) else "dataset"
      paste0(base, "_missing_features_comparison.csv")
    },
    content = function(file) {
      req(rv$result)
      x <- as.data.frame(rv$result$comparison_missing_features)
      readr::write_csv(x, file)
    }
  )
  
  output$dl_nonmissing_csv <- downloadHandler(
    filename = function() {
      base <- if (!is.null(input$file$name)) file_path_sans_ext(input$file$name) else "dataset"
      paste0(base, "_nonmissing_summary.csv")
    },
    content = function(file) {
      req(rv$result)
      x <- as.data.frame(rv$result$summary_nonmissing_features)
      readr::write_csv(x, file)
    }
  )
  
  # --- keep your make_tabset name, but make it a bit sturdier ---
  make_tabset <- function(section, prefix) {
    if (is.null(section) || !length(section)) return(tags$div("No plots yet"))
    
    tabs <- lapply(seq_along(section), function(i) {
      id    <- paste0(prefix, "_", i)
      title <- section[[i]]$variable %||% paste(prefix, i)
      tabPanel(title, plotOutput(id, height = 320))
    })
    
    # renderers for each dynamic plotOutput
    lapply(seq_along(section), function(i) {
      local({
        ii <- i
        id <- paste0(prefix, "_", ii)
        output[[id]] <- renderPlot({
          req(section[[ii]]$combined_plot)
          section[[ii]]$combined_plot
        })
      })
    })
    
    do.call(tabsetPanel, c(list(id = paste0(prefix, "_tabs")), tabs))
  }
  
  
  # --- FIX: use rv$plots instead of plots ---
  output$tabs_missing <- renderUI({
    req(rv$plots)
    make_tabset(rv$plots$comparison_missing_features, "miss")
  })
  
  output$tabs_nonmissing <- renderUI({
    req(rv$plots)
    make_tabset(rv$plots$summary_nonmissing_features, "nonmiss")
  })
  
  output$tabs_missing <- renderUI({
    req(rv$plots)  # FIX: was req(plots)
    validate(need(!is.null(rv$plots$comparison_missing_features), "Run the pipeline to generate plots."))
    make_tabset(rv$plots$comparison_missing_features, "miss")
  })
  
  output$tabs_nonmissing <- renderUI({
    req(rv$plots)  # FIX: was req(plots)
    validate(need(!is.null(rv$plots$summary_nonmissing_features), "Run the pipeline to generate plots."))
    make_tabset(rv$plots$summary_nonmissing_features, "nonmiss")
  })
  
  
  # ---- Session Log ----
  output$log <- renderText({
    paste(rv$log, collapse = "\n")
  })
}

# small infix helper
`%||%` <- function(x, y) if (is.null(x)) y else x

shinyApp(ui, server)
