library(tidyverse)
library(patchwork)

theme_nice <- function(base_size = 12){
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 2),
      plot.subtitle   = element_text(color = "grey30"),
      plot.caption    = element_text(color = "grey40", size = base_size - 2),
      panel.grid.minor= element_blank(),
      legend.position = "top",
      strip.text      = element_text(face = "bold")
    )
}

collapse_levels <- function(x, top_k = 10) {
  x_chr <- as.character(x)
  tab <- sort(table(x_chr, useNA = "no"), decreasing = TRUE)
  keep <- names(head(tab, top_k))
  x_chr[!(x_chr %in% keep) & !is.na(x_chr)] <- "Other"
  factor(x_chr, levels = c(setdiff(keep, NA), "Other"))
}

prep_long_for_var <- function(df_before, df_after, var, cat_top_k = 10) {
  stopifnot(var %in% names(df_before), var %in% names(df_after))
  stopifnot(nrow(df_before) == nrow(df_after))
  
  x_before <- df_before[[var]]
  x_after  <- df_after[[var]]
  was_na   <- is.na(x_before)
  is_num   <- is.numeric(x_after)
  
  if (is_num) {
    tibble(
      variable   = var,
      type       = "numeric",
      source     = c(rep("Before (observed)", sum(!is.na(x_before))),
                     rep("After (imputed)",   length(x_after))),
      value      = c(x_before[!is.na(x_before)], x_after),
      imputedRow = c(rep(FALSE, sum(!is.na(x_before))), was_na)
    )
  } else {
    b_fac <- collapse_levels(x_before, cat_top_k)
    a_fac <- collapse_levels(x_after,  cat_top_k)
    all_lvls <- union(levels(b_fac), levels(a_fac))
    b_fac <- factor(b_fac, levels = all_lvls)
    a_fac <- factor(a_fac, levels = all_lvls)
    
    tibble(
      variable   = var,
      type       = "categorical",
      source     = c(rep("Before (observed)", length(b_fac)),
                     rep("After (imputed)",   length(a_fac))),
      value      = c(b_fac, a_fac),
      imputedRow = c(rep(FALSE, length(b_fac)), was_na)
    )
  }
}

# ===== Plots for comparisons (missing features) =====
plot_numeric_compare <- function(df_long){
  gg_density <- ggplot(df_long, aes(x = value, fill = source)) +
    geom_density(alpha = 0.35) +
    geom_rug(data = subset(df_long, source == "After (imputed)" & imputedRow),
             aes(x = value), sides = "b", alpha = 0.4) +
    labs(title = paste0("Density — ", unique(df_long$variable)),
         subtitle = "(Before vs After)", x = unique(df_long$variable), y = "Density") +
    theme_nice()
  gg_box <- ggplot(df_long, aes(x = source, y = value, fill = source)) +
    geom_boxplot(width = 0.6, outlier.alpha = 0.35) +
    labs(title = "Boxplot", x = NULL, y = unique(df_long$variable)) +
    theme_nice() + theme(legend.position = "none")
  (gg_density | gg_box) + plot_layout(widths = c(2,1))
}

plot_cat_compare <- function(df_long){
  df_plot <- df_long %>%
    filter(!is.na(value)) %>%
    count(source, value, name="n") %>%
    group_by(source) %>% mutate(prop = n/sum(n)) %>% ungroup()
  ggplot(df_plot, aes(x=value, y=prop, fill=source)) +
    geom_col(position=position_dodge(width=0.8), width=0.7) +
    scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
    labs(title=paste0("Category Proportions — ", unique(df_long$variable)),
         subtitle="(Before vs After)", x=unique(df_long$variable), y="Proportion") +
    theme_nice() +
    theme(axis.text.x=element_text(angle=30,hjust=1))
}

# ===== Plots for non-missing features (after only) =====
plot_numeric_single <- function(df_after, var){
  gg_density <- ggplot(df_after, aes(x = .data[[var]])) +
    geom_density(fill="#4DBBD5B2", alpha=0.4) +
    labs(title=paste0("Density — ", var), subtitle="(Non-Missing Feature)",
         x=var, y="Density") + theme_nice()
  gg_box <- ggplot(df_after, aes(y = .data[[var]])) +
    geom_boxplot(fill="#4DBBD5B2", width=0.3, outlier.alpha=0.4) +
    labs(title="Boxplot", x=NULL, y=var) + theme_nice() +
    theme(plot.title=element_text(size=11, face="bold"))
  (gg_density | gg_box) + plot_layout(widths=c(2,1))
}

plot_cat_single <- function(df_after, var, cat_top_k=10){
  f <- collapse_levels(df_after[[var]], cat_top_k)
  df_plot <- tibble(value=f) %>% count(value, name="n") %>% mutate(prop=n/sum(n))
  ggplot(df_plot, aes(x=value, y=prop)) +
    geom_col(fill="#4DBBD5B2", width=0.7) +
    scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
    labs(title=paste0("Category Proportions — ", var),
         subtitle="(Non-Missing Feature)",
         x=var, y="Proportion") +
    theme_nice() +
    theme(axis.text.x=element_text(angle=30,hjust=1))
}

# ===== Main wrapper =====
plot_distributions_comparison <- function(df_before, df_after, cat_top_k=10){
  common <- intersect(names(df_before), names(df_after))
  df_before <- df_before[, common, drop=FALSE]
  df_after  <- df_after[,  common, drop=FALSE]
  
  missing_before <- names(df_before)[colSums(is.na(df_before))>0]
  missing_after  <- names(df_after)[colSums(is.na(df_after))>0]
  missing_features <- union(missing_before, missing_after)
  nonmissing_features <- setdiff(common, missing_features)
  
  # ---- 1️⃣ missing features (Before vs After) ----
  compare_plots <- map(missing_features, function(v){
    df_long <- prep_long_for_var(df_before, df_after, v, cat_top_k)
    if(unique(df_long$type)=="numeric") 
      list(variable=v, combined_plot=plot_numeric_compare(df_long))
    else 
      list(variable=v, combined_plot=plot_cat_compare(df_long))
  })
  
  # ---- 2️⃣ non-missing features (After only) ----
  single_plots <- map(nonmissing_features, function(v){
    if(is.numeric(df_after[[v]]))
      list(variable=v, combined_plot=plot_numeric_single(df_after, v))
    else
      list(variable=v, combined_plot=plot_cat_single(df_after, v, cat_top_k))
  })
  
  list(
    comparison_missing_features = compare_plots,
    summary_nonmissing_features = single_plots
  )
}

# ==========================
# FUNCTION DEFINITION
# ==========================
missing_summary_plot <- function(df) {
  
  # Ensure input is a data frame
  if (!is.data.frame(df)) stop("Input must be a data frame.")
  
  # 1️⃣ Barplot: Number of missing values per feature
  missing_count <- colSums(is.na(df))
  missing_df <- data.frame(Feature = names(missing_count), MissingCount = missing_count)
  
  p1 <- ggplot(missing_df, aes(x = reorder(Feature, -MissingCount), y = MissingCount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = MissingCount), hjust = -0.1, size = 4) +
    coord_flip() +
    labs(title = "Number of Missing Values per Feature",
         x = "Features", y = "Count of Missing Values") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  prop_df <- df %>%
    summarise(across(everything(), ~mean(is.na(.)))) %>%
    pivot_longer(cols = everything(), names_to = "Feature", values_to = "MissingProp") %>%
    mutate(NonMissingProp = 1 - MissingProp) %>%
    pivot_longer(cols = c(NonMissingProp, MissingProp),
                 names_to = "Status", values_to = "Proportion") %>%
    mutate(Status = factor(Status, levels = c("MissingProp", "NonMissingProp")))
  
  p2 <- ggplot(prop_df, aes(x = Feature, y = Proportion, fill = Status)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.9) +
    scale_fill_manual(values = c("MissingProp" = "red", "NonMissingProp" = "gray80"),
                      labels = c("Missing", "Non-Missing")) +
    labs(title = "Proportion of Missing vs Non-Missing per Feature",
         x = "Features", y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  # 3️⃣ Combine plots
  combined_plot <- p1 / p2 + plot_layout(heights = c(1, 1.2))
  
  # Return the final plot
  return(combined_plot)
}

# ==========================
# ✅ Example Usage
# ==========================

summarize_distributions_comparison <- function(before_data, after_data, top_k = 5) {
  
  # Identify missing features
  missing_before <- names(before_data)[colSums(is.na(before_data)) > 0]
  missing_after  <- names(after_data)[colSums(is.na(after_data)) > 0]
  missing_features <- union(missing_before, missing_after)
  
  # Non-missing features
  nonmissing_features <- setdiff(intersect(names(before_data), names(after_data)), missing_features)
  
  # Helper to summarize one dataset
  summarize_data <- function(data, data_name) {
    num_vars <- names(data)[vapply(data, is.numeric, logical(1))]
    cat_vars <- setdiff(names(data), num_vars)
    
    num_tbl <- purrr::map_dfr(num_vars, function(v) {
      x <- data[[v]]
      tibble(
        variable = v, type = "numeric", dataset = data_name,
        n = length(x), n_nonmiss = sum(!is.na(x)),
        mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE),
        min = suppressWarnings(min(x, na.rm = TRUE)),
        max = suppressWarnings(max(x, na.rm = TRUE)),
        skewness = tryCatch(e1071::skewness(x, na.rm = TRUE, type = 2), error = \(e) NA_real_),
        kurtosis = tryCatch(e1071::kurtosis(x, na.rm = TRUE, type = 2), error = \(e) NA_real_)
      )
    })
    
    cat_tbl <- purrr::map_dfr(cat_vars, function(v) {
      x <- as.character(data[[v]])
      lv <- sort(table(x), decreasing = TRUE)
      tibble(
        variable = v, type = "categorical", dataset = data_name,
        n = length(x), n_nonmiss = sum(!is.na(x)),
        levels = paste0(names(head(lv, top_k)), collapse = " | "),
        counts = paste0(as.integer(head(lv, top_k)), collapse = " | ")
      )
    })
    
    bind_rows(num_tbl, cat_tbl)
  }
  
  # ---- Part 1: Comparison for Missing Features ----
  before_summary <- summarize_data(before_data[, missing_features, drop = FALSE], "Before")
  after_summary  <- summarize_data(after_data[, missing_features, drop = FALSE], "After")
  
  comparison_tbl <- full_join(before_summary, after_summary,
                              by = c("variable", "type"),
                              suffix = c("_before", "_after"))
  
  # ---- Part 2: Summary for Non-Missing Features ----
  nonmissing_summary <- summarize_data(after_data[, nonmissing_features, drop = FALSE], "After (NonMissing)")
  
  list(
    comparison_missing_features = comparison_tbl,
    summary_nonmissing_features = nonmissing_summary
  )
}

correlation_compare_plots <- function(before_data,
                                      after_data,
                                      vars = NULL,
                                      method = c("pearson","spearman","kendall"),
                                      min_var = 1e-12,
                                      reorder = TRUE,
                                      title_prefix = "Correlation") {
  method <- match.arg(method)
  
  stopifnot(is.data.frame(before_data), is.data.frame(after_data))
  common <- intersect(names(before_data), names(after_data))
  if (!length(common)) stop("No common columns between before and after.")
  
  # numeric-only, shared columns
  num_in_b <- names(which(vapply(before_data[common], is.numeric, logical(1))))
  num_in_a <- names(which(vapply(after_data[common],  is.numeric, logical(1))))
  num_vars <- intersect(num_in_b, num_in_a)
  
  if (!is.null(vars)) {
    # user-provided subset: keep only those that are numeric in both
    vars <- intersect(vars, num_vars)
    num_vars <- vars
  }
  
  if (length(num_vars) < 2) {
    stop("Fewer than 2 continuous (numeric) features in common; cannot compute correlations.")
  }
  
  B <- before_data[, num_vars, drop = FALSE]
  A <- after_data[,  num_vars, drop = FALSE]
  
  # drop near-constant columns (variance below threshold) *in either* dataset
  varB <- vapply(B, function(x) var(x, na.rm = TRUE), numeric(1))
  varA <- vapply(A, function(x) var(x, na.rm = TRUE), numeric(1))
  keep <- names(which(varB > min_var & varA > min_var))
  if (length(keep) < 2) stop("After removing near-constant columns, fewer than 2 numeric variables remain.")
  B <- B[, keep, drop = FALSE]
  A <- A[, keep, drop = FALSE]
  
  # compute correlations
  cor_before <- suppressWarnings(stats::cor(B, use = "pairwise.complete.obs", method = method))
  cor_after  <- suppressWarnings(stats::cor(A, use = "pairwise.complete.obs", method = method))
  # align (defensive)
  cor_after  <- cor_after[rownames(cor_before), colnames(cor_before)]
  cor_delta  <- cor_after - cor_before
  
  # optional reordering by hierarchical clustering on 1 - |r|
  reorder_idx <- function(C) {
    D <- as.dist(1 - abs(C))
    hc <- hclust(D, method = "average")
    hc$order
  }
  if (reorder) {
    ord <- reorder_idx(cor_after)  # order by AFTER (often cleaner)
    cor_before <- cor_before[ord, ord, drop = FALSE]
    cor_after  <- cor_after [ord, ord, drop = FALSE]
    cor_delta  <- cor_delta [ord, ord, drop = FALSE]
  }
  
  # ---- plotting helper (ggplot heatmap) ----
  mat_to_long <- function(M) {
    tibble::tibble(
      row = rep(rownames(M), times = ncol(M)),
      col = rep(colnames(M), each  = nrow(M)),
      r   = as.vector(M)
    )
  }
  
  lim1 <- 1  # for correlations
  limD <- max(0.1, max(abs(c(cor_delta)), na.rm = TRUE))  # symmetric limits for delta
  
  plot_mat <- function(M, title, lim = 1) {
    df <- mat_to_long(M)
    # keep full matrix (diagonal included), it's often useful for scan
    ggplot(df, aes(x = col, y = row, fill = r)) +
      geom_tile(color = "white", size = 0.2) +
      scale_fill_gradient2(limits = c(-lim, lim),
                           low = "#4575b4", mid = "white", high = "#d73027",
                           name = "r") +
      coord_fixed() +
      labs(title = title, x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold")
      )
  }
  
  p_before <- plot_mat(cor_before, sprintf("%s — Before (%s)", title_prefix, method), lim = lim1)
  p_after  <- plot_mat(cor_after,  sprintf("%s — After  (%s)", title_prefix, method), lim = lim1)
  p_delta  <- plot_mat(cor_delta,  sprintf("Δ After − Before (%s)", method),           lim = limD)
  
  combined <- (p_before | p_after | p_delta) + patchwork::plot_layout(widths = c(1,1,1))
  
  list(
    mats = list(before = cor_before, after = cor_after, delta = cor_delta),
    plots = list(before = p_before, after = p_after, delta = p_delta, combined = combined),
    vars_used = colnames(cor_after)
  )
}

# missing_summary_plot(df)

# imputed_data <- out$best_imputed[[5]]
# 
# result <- summarize_distributions_comparison(before_data = df, after_data = imputed_data)
# 
# # View results
# data.frame(result$comparison_missing_features)   # Comparison table (before vs after)
# data.frame(result$summary_nonmissing_features)   # Summary table for non-missing features
# 
# 
# 
# plots <- plot_distributions_comparison(df, imputed_data, cat_top_k = 8)
# 
# plots
# 
# cc <- correlation_compare_plots(df, imputed_data, method = "pearson", reorder = TRUE)
# 
# # Show all three panels: Before | After | Δ
# print(cc$plots$combined)

