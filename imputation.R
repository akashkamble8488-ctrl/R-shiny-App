# ============ Unified Iterative Imputer (6 methods) ============
# Methods:
# 1) method = "pmm"          -> Predictive Mean Matching (numeric targets)
# 2) method = "linreg"       -> OLS (numeric) / logistic (binary); >2 levels -> CART fallback
# 3) method = "linreg_boot"  -> Bootstrap OLS / logistic; >2 levels -> CART fallback
# 4) method = "bayes"        -> Conjugate Bayes linear reg (numeric only); categorical -> CART fallback
# 5) method = "cart"         -> rpart
# 6) method = "rf"           -> randomForest

impute_iterative <- function(data,
                             method = c("rf","cart","pmm","linreg","linreg_boot","bayes"),
                             m = 1,
                             max_iter = 5,
                             tol = 1e-3,
                             seed = 123,
                             ntree = 300,         # for RF
                             k_pmm = 5,           # for PMM donors
                             verbose = TRUE) {
  method <- match.arg(method)
  stopifnot(is.data.frame(data))
  set.seed(seed)
  
  # --- Dependencies (conditionally) ---
  if (method == "rf")  { if (!requireNamespace("randomForest", quietly=TRUE)) stop("Need randomForest") }
  if (method == "cart"){ if (!requireNamespace("rpart", quietly=TRUE))        stop("Need rpart") }
  if (method == "bayes"){ if (!requireNamespace("MASS", quietly=TRUE))        stop("Need MASS") }
  
  # --- Helpers ---
  get_mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }
  is_binary_factor <- function(x) is.factor(x) && length(levels(x)) == 2
  
  # Build a modeling dataset for a given target: return list(train_x, train_y, new_x)
  build_xy <- function(d, target, na_idx, obs_idx, predictors){
    train_x <- d[obs_idx, predictors, drop = FALSE]
    train_y <- d[obs_idx, target, drop = FALSE][[1]]
    new_x   <- d[na_idx, predictors, drop = FALSE]
    
    # Clean residual NAs in predictors (use current d stats)
    for (p in predictors) {
      if (anyNA(train_x[[p]])) {
        if (is.numeric(train_x[[p]])) {
          train_x[[p]][is.na(train_x[[p]])] <- mean(d[[p]], na.rm = TRUE)
        } else {
          train_x[[p]][is.na(train_x[[p]])] <- get_mode(d[[p]])
        }
      }
      if (anyNA(new_x[[p]])) {
        if (is.numeric(new_x[[p]])) {
          new_x[[p]][is.na(new_x[[p]])] <- mean(d[[p]], na.rm = TRUE)
        } else {
          new_x[[p]][is.na(new_x[[p]])] <- get_mode(d[[p]])
        }
      }
    }
    list(train_x=train_x, train_y=train_y, new_x=new_x)
  }
  
  # Fit/predict for each method ------------------------------------------------
  fit_predict <- function(method, d, target, predictors, na_idx, obs_idx){
    if (length(na_idx) == 0 || length(obs_idx) < 5) return(NULL)
    
    xy <- build_xy(d, target, na_idx, obs_idx, predictors)
    train_x <- xy$train_x; train_y <- xy$train_y; new_x <- xy$new_x
    y_is_num <- is.numeric(d[[target]])
    
    # --- PMM (numeric only) ---
    if (method == "pmm") {
      if (!y_is_num) return(list(pred = d[[target]][na_idx]))  # fallback: no change for non-numeric
      # Linear model on observed
      df_train <- data.frame(.y = train_y, train_x, check.names = FALSE)
      fit <- tryCatch(stats::lm(.y ~ ., data = df_train), error=function(e) NULL)
      if (is.null(fit)) return(NULL)
      
      preds_obs <- as.numeric(stats::predict(fit, newdata = train_x))
      preds_mis <- as.numeric(stats::predict(fit, newdata = new_x))
      y_obs     <- as.numeric(train_y)
      
      imputed <- numeric(length(na_idx))
      for (i in seq_along(na_idx)) {
        dists <- abs(preds_obs - preds_mis[i])
        donors <- order(dists)[seq_len(min(k_pmm, length(dists)))]
        imputed[i] <- sample(y_obs[donors], 1)
      }
      return(list(pred = imputed))
    }
    
    # --- Linear regression / bootstrap / Bayes ---
    if (method %in% c("linreg","linreg_boot","bayes")) {
      # CATEGORICAL branch:
      if (!y_is_num) {
        # If binary: logistic (glm), else CART fallback
        if (is_binary_factor(d[[target]])) {
          lvl <- levels(d[[target]])
          df_train <- data.frame(.y = factor(train_y, levels=lvl), train_x, check.names = FALSE)
          
          if (method == "linreg_boot") {
            boot_rows <- sample(seq_len(nrow(df_train)), nrow(df_train), replace = TRUE)
            df_train <- df_train[boot_rows, , drop=FALSE]
          }
          fit <- tryCatch(stats::glm(.y ~ ., data = df_train, family = binomial()), error=function(e) NULL)
          if (is.null(fit)) return(NULL)
          p_hat <- stats::predict(fit, newdata = new_x, type = "response")
          # sample class using Bernoulli
          draw <- rbinom(length(p_hat), 1, p_hat)
          pred <- factor(ifelse(draw==1, lvl[2], lvl[1]), levels=lvl)
          return(list(pred = pred))
        } else {
          # Multiclass fallback -> CART
          method_cart <- "cart"
          # Fall through to CART below via recursion
          return(fit_predict(method_cart, d, target, predictors, na_idx, obs_idx))
        }
      }
      
      # NUMERIC branch:
      if (method == "linreg") {
        df_train <- data.frame(.y = train_y, train_x, check.names = FALSE)
        fit <- tryCatch(stats::lm(.y ~ ., data = df_train), error=function(e) NULL)
        if (is.null(fit)) return(NULL)
        preds <- as.numeric(stats::predict(fit, newdata = new_x))
        # Add residual noise to avoid underestimating variance
        sig  <- summary(fit)$sigma
        return(list(pred = preds + rnorm(length(preds), 0, sig)))
      }
      
      if (method == "linreg_boot") {
        df_train <- data.frame(.y = train_y, train_x, check.names = FALSE)
        boot_rows <- sample(seq_len(nrow(df_train)), nrow(df_train), replace = TRUE)
        dfb <- df_train[boot_rows, , drop=FALSE]
        fit <- tryCatch(stats::lm(.y ~ ., data = dfb), error=function(e) NULL)
        if (is.null(fit)) return(NULL)
        preds <- as.numeric(stats::predict(fit, newdata = new_x))
        sig  <- summary(fit)$sigma
        return(list(pred = preds + rnorm(length(preds), 0, sig)))
      }
      
      if (method == "bayes") {
        # Conjugate Normal-Inverse-Gamma posterior draw
        # y ~ N(Xb, s2 I), prior vague; beta|s2 ~ N(bhat, s2 (X'X)^-1), s2 ~ scaled-Inv-chi^2
        X <- model.matrix(~ ., data = train_x)
        y <- as.numeric(train_y)
        # OLS
        bhat <- tryCatch(solve(t(X) %*% X, t(X) %*% y), error=function(e) NULL)
        if (is.null(bhat)) {
          # try pseudo-inverse
          XtX <- t(X) %*% X
          bhat <- MASS::ginv(XtX) %*% t(X) %*% y
        }
        resid <- y - as.numeric(X %*% bhat)
        n <- length(y); p <- ncol(X)
        df <- max(1, n - p)
        sigma2_hat <- sum(resid^2)/df
        XtX_inv <- tryCatch(solve(t(X) %*% X), error=function(e) MASS::ginv(t(X) %*% X))
        
        sigma2_draw <- sigma2_hat * df / stats::rchisq(1, df)
        beta_draw   <- as.numeric(MASS::mvrnorm(1, mu = as.numeric(bhat), Sigma = sigma2_draw * XtX_inv))
        Xm <- model.matrix(~ ., data = new_x)
        mu <- as.numeric(Xm %*% beta_draw)
        pred <- mu + stats::rnorm(length(mu), 0, sqrt(sigma2_draw))
        return(list(pred = pred))
      }
    }
    
    # --- CART ---
    if (method == "cart") {
      meth <- if (y_is_num) "anova" else "class"
      df_train <- data.frame(.y = train_y, train_x, check.names = FALSE)
      fit <- tryCatch(rpart::rpart(.y ~ ., data = df_train, method = meth),
                      error=function(e) NULL)
      if (is.null(fit)) return(NULL)
      if (y_is_num) {
        pr <- as.numeric(stats::predict(fit, newdata = new_x))
        # tiny noise to reduce ties/underdispersion
        pr <- pr + rnorm(length(pr), 0, sd(pr, na.rm=TRUE) * 0.01)
        return(list(pred = pr))
      } else {
        pr <- stats::predict(fit, newdata = new_x, type = "prob")
        if (is.matrix(pr)) {
          levs <- colnames(pr)
          draw <- apply(pr, 1, function(p) sample(levs, 1, prob = p))
          return(list(pred = factor(draw, levels = levels(d[[target]]))))
        } else {
          # type="class" fallback
          prc <- stats::predict(fit, newdata = new_x, type = "class")
          return(list(pred = prc))
        }
      }
    }
    
    # --- RF ---
    if (method == "rf") {
      df_train <- data.frame(train_x, .y = train_y, check.names = FALSE)
      form <- stats::reformulate(termlabels = colnames(train_x), response = ".y")
      fit <- tryCatch(
        randomForest::randomForest(form, data = df_train, ntree = ntree, na.action = na.omit),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      pr <- stats::predict(fit, newdata = new_x)
      if (!is.numeric(d[[target]])) pr <- factor(pr, levels = levels(d[[target]]))
      return(list(pred = pr))
    }
    
    NULL
  }
  
  # --- Core worker to produce ONE imputed dataset (single chain) ---------------
  run_one <- function() {
    d <- data
    
    # Convert character -> factor (safer for tree/RF/cart)
    char_cols <- vapply(d, is.character, logical(1))
    d[char_cols] <- lapply(d[char_cols], function(x) factor(x, exclude = NULL))
    
    # Map original NA positions; choose variables with missingness
    na_map <- lapply(d, function(col) which(is.na(col)))
    vars_with_na <- names(na_map)[vapply(na_map, length, integer(1)) > 0]
    if (length(vars_with_na) == 0) {
      if (verbose) message("No missing values.")
      return(list(data = d, history = data.frame(iter=0,total_change=0)))
    }
    
    # Initial simple imputation (mean/mode)
    for (v in vars_with_na) {
      idx <- na_map[[v]]
      if (is.numeric(d[[v]])) {
        d[[v]][idx] <- mean(d[[v]], na.rm = TRUE)
      } else {
        mv <- get_mode(d[[v]])
        d[[v]][idx] <- mv
        d[[v]] <- droplevels(d[[v]])
      }
    }
    if (verbose) message("Initial imputation for: ", paste(vars_with_na, collapse=", "))
    
    history <- data.frame(iter=integer(0), total_change=numeric(0))
    
    for (iter in seq_len(max_iter)) {
      if (verbose) cat("\nIteration", iter, "\n")
      old <- d
      
      for (v in vars_with_na) {
        na_idx  <- na_map[[v]]
        if (length(na_idx) == 0) next
        obs_idx <- setdiff(seq_len(nrow(d)), na_idx)
        predictors <- setdiff(names(d), v)
        
        fp <- fit_predict(method, d, v, predictors, na_idx, obs_idx)
        if (is.null(fp)) next
        d[[v]][na_idx] <- fp$pred
        if (!is.numeric(d[[v]]) && !is.factor(d[[v]])) d[[v]] <- factor(d[[v]])
      }
      
      # Convergence metric: RMSE for numeric (at NA rows), fraction changed for categorical
      total_change <- 0
      for (v in vars_with_na) {
        na_idx <- na_map[[v]]
        if (length(na_idx) == 0) next
        if (is.numeric(d[[v]])) {
          delta <- old[[v]][na_idx] - d[[v]][na_idx]
          total_change <- total_change + sqrt(mean(delta^2, na.rm = TRUE))
        } else {
          changed <- sum(old[[v]][na_idx] != d[[v]][na_idx], na.rm = TRUE) / length(na_idx)
          total_change <- total_change + changed
        }
      }
      history <- rbind(history, data.frame(iter=iter, total_change=total_change))
      if (verbose) cat("   ↳ change:", signif(total_change, 6), "\n")
      if (total_change < tol) {
        if (verbose) cat("Converged at iteration", iter, "\n")
        break
      }
    }
    
    list(data=d, history=history)
  }
  
  # --- Multiple imputations (m datasets) --------------------------------------
  if (m == 1) {
    out <- run_one()
    return(out)  # list(data=..., history=...)
  } else {
    res <- vector("list", m)
    hist <- vector("list", m)
    for (i in seq_len(m)) {
      if (verbose) cat("\n=== Chain", i, "of", m, "===\n")
      # jitter seed across chains
      set.seed(seed + 1000*i)
      ri <- run_one()
      res[[i]] <- ri$data
      hist[[i]] <- ri$history
    }
    return(list(data_list = res, histories = hist))
  }
}

rubin_pool <- function(q_vec, se_vec) {
  stopifnot(length(q_vec) == length(se_vec), length(q_vec) >= 1)
  ok <- is.finite(q_vec) & is.finite(se_vec)
  if (!all(ok)) {
    q_vec <- q_vec[ok]; se_vec <- se_vec[ok]
  }
  m <- length(q_vec)
  if (m == 0) stop("No valid estimates to pool.")
  qbar <- mean(q_vec)
  ubar <- mean(se_vec^2)
  
  b <- if (m == 1) 0 else stats::var(q_vec)
  tvar <- ubar + (1 + 1/m) * b
  r <- if (ubar == 0) Inf else (1 + 1/m) * b / ubar
  df <- if (is.infinite(r)) (m - 1) else (m - 1) * (1 + 1/r)^2
  
  list(
    qbar = qbar,
    se   = sqrt(tvar),
    df   = df,
    ubar = ubar,
    b    = b,
    t    = tvar
  )
}

run_all_imputations <- function(data, 
                                methods = c("pmm","linreg","linreg_boot","bayes","cart","rf"),
                                m = 5, max_iter = 6, tol = 1e-3, seed = 123, ntree = 500, k_pmm = 5, verbose = TRUE, ...) {
  
  results <- vector("list", length(methods)); names(results) <- methods
  for (i in seq_along(methods)) {
    meth <- methods[i]
    if (verbose) cat("\n==============================\nMethod:", toupper(meth), "\n==============================\n")
    set.seed(seed + i * 10000L)  # de-correlate per method
    results[[i]] <- tryCatch(
      impute_iterative(
        data, method = meth, m = m, max_iter = max_iter, tol = tol,
        seed = seed + i * 1000L, ntree = ntree, k_pmm = k_pmm,
        verbose = verbose, ...
      ),
      error = function(e) { if (verbose) message("  ✗ ", meth, ": ", e$message); NULL }
    )
  }
  results
}


pool_term <- function(res, formula, term, fit_fun = stats::lm) {
  data_list <- if (!is.null(res$data_list)) res$data_list else if (!is.null(res$data)) list(res$data) else stop("Unexpected structure.")
  q <- se <- rep(NA_real_, length(data_list))
  
  for (i in seq_along(data_list)) {
    d <- data_list[[i]]
    fit <- tryCatch(fit_fun(formula, data = d), error = function(e) NULL)
    if (is.null(fit)) next
    co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(co) || !(term %in% rownames(co))) next
    q[i]  <- co[term, "Estimate"]
    se[i] <- co[term, "Std. Error"]
  }
  
  if (all(is.na(q))) stop("No valid estimates for term '", term, "'.")
  rubin_pool(q[is.finite(q) & is.finite(se)], se[is.finite(q) & is.finite(se)])
}


## ============================================================
## Auto-impute & select best method (works for mixed data)
## - Detects missing columns
## - Evaluates each method via repeated hold-out masking
## - Ranks and returns best result + a leaderboard
## Requires: impute_iterative(), run_all_imputations()
## ============================================================

auto_impute_select <- function(data,
                               methods    = c("pmm","linreg","linreg_boot","bayes","cart","rf"),
                               m          = 5,
                               max_iter   = 6,
                               tol        = 1e-3,
                               seed       = 123,
                               ntree      = 500,
                               k_pmm      = 5,
                               eval_frac  = 0.10,   # % of observed values to mask per column with NA
                               eval_reps  = 3,      # repeats for stable estimates
                               verbose    = TRUE) {
  stopifnot(is.data.frame(data))
  set.seed(seed)
  
  # ---------- 0) Identify missing columns ----------
  na_map <- lapply(data, function(x) which(is.na(x)))
  miss_cols <- names(na_map)[vapply(na_map, length, integer(1)) > 0]
  if (length(miss_cols) == 0) {
    if (verbose) message("No missing values detected. Returning original data.")
    return(list(
      original         = data,
      missing_map      = na_map,
      best_method      = NA_character_,
      best_imputed     = data,
      leaderboard      = data.frame(),
      method_results   = list(),
      histories        = list()
    ))
  }
  if (verbose) message("Missing columns: ", paste(miss_cols, collapse=", "))
  
  # Helpers
  is_cat <- function(x) is.factor(x) || is.character(x)
  is_num <- function(x) is.numeric(x)
  to_factor <- function(d) {
    d2 <- d
    char_cols <- vapply(d2, is.character, logical(1))
    d2[char_cols] <- lapply(d2[char_cols], function(x) factor(x, exclude=NULL))
    d2
  }
  
  # ---------- 1) Evaluation-by-masking metrics ----------
  metric_one <- function(col, truth, pred) {
    if (is_num(truth)) {
      rmse <- sqrt(mean((truth - pred)^2, na.rm=TRUE))
      c(num_rmse = rmse, cat_acc = NA_real_)
    } else {
      # ensure factors comparable
      acc <- mean(as.character(truth) == as.character(pred))
      c(num_rmse = NA_real_, cat_acc = acc)
    }
  }
  
  eval_method_once <- function(method) {
    # One evaluation repetition across all cols with missingness:
    d_eval <- data
    masks <- list()
    for (v in miss_cols) {
      obs_idx <- which(!is.na(d_eval[[v]]))
      if (length(obs_idx) == 0) next
      k <- max(1, floor(eval_frac * length(obs_idx)))
      masks[[v]] <- sample(obs_idx, k)
      d_eval[[v]][masks[[v]]] <- NA
    }
    # Impute once (m=1) for speed in eval
    res <- tryCatch(
      impute_iterative(
        data = to_factor(d_eval),
        method = method, m = 1, max_iter = max_iter, tol = tol,
        seed = seed + sample.int(1e6, 1),
        ntree = ntree, k_pmm = k_pmm, verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(res)) return(NULL)
    di <- res$data
    # Collect metrics across masked cells
    mrows <- lapply(miss_cols, function(v) {
      idx <- masks[[v]]; if (is.null(idx) || length(idx)==0) return(NULL)
      truth <- data[[v]][idx]
      pred  <- di[[v]][idx]
      as.list(metric_one(v, truth, pred))
    })
    # Aggregate
    num_rmse <- mean(sapply(mrows, function(z) z$num_rmse), na.rm=TRUE)
    cat_acc  <- mean(sapply(mrows, function(z) z$cat_acc),  na.rm=TRUE)
    list(num_rmse = num_rmse, cat_acc = cat_acc)
  }
  
  eval_method <- function(method) {
    # Repeat to stabilize
    scores <- replicate(eval_reps, eval_method_once(method), simplify = FALSE)
    scores <- scores[!vapply(scores, is.null, logical(1))]
    if (!length(scores)) return(c(num_rmse=Inf, cat_acc=0))
    c(
      num_rmse = mean(vapply(scores, `[[`, numeric(1), "num_rmse"), na.rm=TRUE),
      cat_acc  = mean(vapply(scores, `[[`, numeric(1), "cat_acc"),  na.rm=TRUE)
    )
  }
  
  # ---------- 2) Full imputation runs (for outputs) ----------
  # Also track runtime and convergence
  runtimes <- numeric(length(methods))
  results  <- vector("list", length(methods)); names(results) <- methods
  histories <- vector("list", length(methods)); names(histories) <- methods
  
  if (verbose) cat("\n=== Running full imputations for all methods ===\n")
  for (i in seq_along(methods)) {
    meth <- methods[i]
    if (verbose) cat("\n----", toupper(meth), "----\n")
    t0 <- proc.time()[3]
    results[[i]] <- tryCatch(
      impute_iterative(
        data = to_factor(data),
        method = meth, m = m, max_iter = max_iter, tol = tol,
        seed = seed + i*1000L, ntree = ntree, k_pmm = k_pmm, verbose = verbose
      ),
      error = function(e) { if (verbose) message("  ✗ ", e$message); NULL }
    )
    runtimes[i] <- proc.time()[3] - t0
    histories[[i]] <- tryCatch(results[[i]]$history %||% results[[i]]$histories, silent = TRUE)
  }
  
  # ---------- 3) Objective evaluation & ranking ----------
  if (verbose) cat("\n=== Evaluating methods by masking ===\n")
  
  # Build a named list of per-method evals; fill failures with NA
  eval_list <- setNames(vector("list", length(methods)), methods)
  for (mm in methods) {
    sc <- tryCatch(eval_method(mm), error = function(e) NULL)
    if (is.null(sc) || !length(sc)) {
      eval_list[[mm]] <- list(num_rmse = NA_real_, cat_acc = NA_real_)
    } else {
      # ensure both fields exist
      nr <- if (!is.null(sc[["num_rmse"]])) sc[["num_rmse"]] else NA_real_
      ca <- if (!is.null(sc[["cat_acc"]]))  sc[["cat_acc"]]  else NA_real_
      eval_list[[mm]] <- list(num_rmse = nr, cat_acc = ca)
    }
  }
  
  # Bind to a data frame with a 'method' column; no rownames assignment
  eval_df <- do.call(
    rbind,
    lapply(names(eval_list), function(nm) {
      data.frame(method = nm,
                 num_rmse = eval_list[[nm]]$num_rmse,
                 cat_acc  = eval_list[[nm]]$cat_acc,
                 stringsAsFactors = FALSE)
    })
  )
  
  # Normalize metrics (lower RMSE better; higher ACC better)
  num_vals <- eval_df$num_rmse
  cat_vals <- eval_df$cat_acc
  
  zrmse <- if (all(is.na(num_vals))) rep(0, nrow(eval_df)) else {
    rng <- range(num_vals, na.rm = TRUE)
    (num_vals - rng[1]) / (diff(rng) + 1e-12)
  }
  
  z1macc <- if (all(is.na(cat_vals))) rep(0, nrow(eval_df)) else {
    v <- 1 - cat_vals
    rng <- range(v, na.rm = TRUE)
    (v - rng[1]) / (diff(rng) + 1e-12)
  }
  
  # Dataset-type weights
  has_num <- any(!is.na(num_vals))
  has_cat <- any(!is.na(cat_vals))
  w_num <- if (has_num) 0.5 else 0
  w_cat <- if (has_cat) 0.5 else 0
  if (w_num + w_cat == 0) { w_num <- 0.5; w_cat <- 0.5 }  # fallback
  
  total_score <- w_num * zrmse + w_cat * z1macc
  
  # Convergence proxy: last total_change
  last_change <- sapply(methods, function(mm) {
    r <- results[[mm]]
    if (is.null(r)) return(Inf)
    if (!is.null(r$history) && nrow(r$history)) {
      tail(r$history$total_change, 1)
    } else if (!is.null(r$histories) && length(r$histories) && nrow(r$histories[[1]])) {
      tail(r$histories[[1]]$total_change, 1)
    } else NA_real_
  })
  
  # Runtime vector already aligned by 'methods'
  leaderboard <- data.frame(
    method       = eval_df$method,
    num_rmse     = eval_df$num_rmse,
    cat_acc      = eval_df$cat_acc,
    total_score  = as.numeric(total_score),
    runtime_sec  = round(runtimes, 3),
    last_change  = last_change,
    stringsAsFactors = FALSE
  )
  leaderboard <- leaderboard[order(leaderboard$total_score, leaderboard$runtime_sec), ]
  
  
  # ---------- 4) Pick best & materialize best imputed dataset ----------
  best_method <- leaderboard$method[1]
  best_res <- results[[best_method]]
  # Choose a representative dataset:
  best_imputed <- if (!is.null(best_res$data)) best_res$data else best_res$data_list[[1]]
  
  # ---------- 5) Return everything ----------
  list(
    original        = data,
    missing_map     = na_map,
    best_method     = best_method,
    best_imputed    = best_imputed,
    leaderboard     = leaderboard,
    method_results  = results,
    histories       = histories
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

auto_clean_mixed_data <- function(df) {
  # --- 1️⃣ Detect ID-like columns ---
  # Logic: likely ID if all values are unique OR name includes id-like pattern
  id_patterns <- names(df)[sapply(df, is.character)]
  
  id_cols <- names(which(sapply(df[,id_patterns], function(x) length(table(x)))/nrow(df)>0.9))
  
  # --- 2️⃣ Remove detected ID columns ---
  df <- df[, !(names(df) %in% id_cols), drop = FALSE]
  
  # --- 3️⃣ Convert remaining character columns to factors ---
  for (col in names(df)) {
    if (is.character(df[[col]])) {
      df[[col]] <- as.factor(df[[col]])
    }
  }
  
  # --- 4️⃣ Ensure numeric columns stay numeric ---
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) return(x)
    suppressWarnings(as.numeric(x))
  })
  
  # --- 5️⃣ Return cleaned data ---
  return(df)
}

# data <- read.csv("C:/Users/tejas/Downloads/healthcare_dataset_with_missing.csv")
# 
# 
# df=auto_clean_mixed_data(data)
# 
# set.seed(123)
# out <- auto_impute_select(
#   data       = df,
#   methods    = c("pmm","linreg","linreg_boot","bayes","cart","rf"),
#   m          = 5,
#   max_iter   = 10,
#   eval_frac  = 0.1,
#   eval_reps  = 3,
#   verbose    = TRUE
# )
# 
# # Outputs
# out$best_method          
# head(out$best_imputed)   # your chosen, fully imputed dataset
# out$leaderboard          # per-method metrics table (RMSE/Acc/score/runtime/convergence)

