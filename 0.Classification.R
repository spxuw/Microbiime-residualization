# ============================================
# Per-BioProject modeling (BINARY Phenotype only)
# Targets:
#   (A) Phenotype (binary)  -> Null: Age-only, Full: Taxa+Age (RF + XGB)
#   (B) AgeGroup (binary)   -> RF + XGB (optional; kept)
# CV: repeated stratified K-fold
# Metrics: AUROC; Key output: Delta AUC = AUC(Full) - AUC(Null) for Phenotype
# Importance: averaged across CV fits (NO full-data model)
# ============================================

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
  library(ranger)
  library(xgboost)
})

setwd("/Users/xuwenwang/Dropbox/Projects/Deage/code")

# ---------- helpers ----------

make_age_group_median <- function(age) {
  thr <- stats::median(age, na.rm = TRUE)
  factor(ifelse(age >= thr, "Older", "Younger"), levels = c("Younger", "Older"))
}

make_folds <- function(y, k = 5, seed = 1) {
  set.seed(seed)
  y <- as.factor(y)
  idx <- seq_along(y)
  folds <- vector("list", k)
  
  for (lv in levels(y)) {
    ii <- idx[y == lv]
    ii <- sample(ii)
    split_ii <- split(ii, rep_len(1:k, length(ii)))
    for (j in 1:k) folds[[j]] <- c(folds[[j]], split_ii[[j]])
  }
  folds
}

# Binary AUROC (expects prob vector for "positive" class = levels(y)[2])
auc_binary <- function(y_true, p_pos) {
  y_true <- factor(y_true)
  if (nlevels(y_true) != 2) return(NA_real_)
  if (length(unique(y_true)) < 2) return(NA_real_)
  roc_obj <- pROC::roc(y_true, as.numeric(p_pos), quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}

drop_allzero_features <- function(df, meta_cols = c("Age", "Phenotype", "BioProject")) {
  feat_cols <- setdiff(colnames(df), meta_cols)
  X <- df[, feat_cols, drop = FALSE]
  X_mat <- as.matrix(X)
  storage.mode(X_mat) <- "numeric"
  keep <- colSums(X_mat != 0, na.rm = TRUE) > 0
  df[, c(feat_cols[keep], meta_cols), drop = FALSE]
}

add_named <- function(acc, x) {
  if (is.null(x) || length(x) == 0) return(acc)
  if (is.null(acc)) return(x)
  alln <- union(names(acc), names(x))
  acc2 <- setNames(rep(0, length(alln)), alln)
  acc2[names(acc)] <- acc2[names(acc)] + acc
  acc2[names(x)]   <- acc2[names(x)] + x
  acc2
}

rf_imp_vec <- function(rf_model) {
  imp <- ranger::importance(rf_model)
  as.numeric(imp) |> setNames(names(imp))
}

xgb_imp_vec <- function(xgb_model, feature_names) {
  imp <- xgboost::xgb.importance(model = xgb_model, feature_names = feature_names)
  if (is.null(imp) || nrow(imp) == 0) return(NULL)
  v <- imp$Gain
  names(v) <- imp$Feature
  v
}

# ---------- core runner for one BioProject (BINARY Phenotype only) ----------
cv_mean_delta_auc <- function(y_pheno, X_taxa, age_vec, k, repeats, seed) {
  # y_pheno must be factor with exactly 2 levels
  stopifnot(nlevels(y_pheno) == 2)
  
  delta_rf <- c()
  delta_xgb <- c()
  
  for (r in seq_len(repeats)) {
    folds <- make_folds(y_pheno, k = k, seed = seed + 1000*r)
    
    for (i in seq_len(k)) {
      te <- folds[[i]]
      tr <- setdiff(seq_along(y_pheno), te)
      
      ytr <- droplevels(y_pheno[tr])
      yte <- droplevels(y_pheno[te])
      if (nlevels(ytr) != 2 || nlevels(yte) != 2) next
      if (any(table(ytr) < 2)) next
      
      Xtr_taxa <- X_taxa[tr, , drop = FALSE]
      Xte_taxa <- X_taxa[te, , drop = FALSE]
      age_tr <- age_vec[tr]
      age_te <- age_vec[te]
      
      # null (age-only)
      lev <- levels(ytr)
      ytr_bin <- ifelse(ytr == lev[2], 1, 0)
      fit_null <- glm(ytr_bin ~ age_tr, family = binomial())
      p_null <- as.numeric(predict(fit_null, newdata = data.frame(age_tr = age_te), type = "response"))
      auc0 <- auc_binary(yte, p_null)
      if (!is.finite(auc0)) next
      
      # full RF (taxa+age)
      Xtr_full <- cbind(Xtr_taxa, Age = age_tr)
      Xte_full <- cbind(Xte_taxa, Age = age_te)
      
      rf <- ranger::ranger(
        dependent.variable.name = "y",
        data = data.frame(y = ytr, Xtr_full, check.names = FALSE),
        probability = TRUE,
        num.trees = 500,
        seed = seed + 10*r + i
      )
      prf <- as.matrix(predict(rf, data = data.frame(Xte_full, check.names = FALSE))$predictions)
      p_rf <- if (lev[2] %in% colnames(prf)) prf[, lev[2]] else prf[, 2]
      auc_rf <- auc_binary(yte, p_rf)
      
      # full XGB (taxa+age)
      dtr <- xgb.DMatrix(data = as.matrix(Xtr_full), label = ytr_bin)
      dte <- xgb.DMatrix(data = as.matrix(Xte_full))
      
      xgb_fit <- xgb.train(
        params = list(
          objective = "binary:logistic",
          eval_metric = "auc",
          max_depth = 6,
          eta = 0.05,
          subsample = 0.8,
          colsample_bytree = 0.8
        ),
        data = dtr,
        nrounds = 300,
        verbose = 0
      )
      p_xgb <- as.numeric(predict(xgb_fit, dte))
      auc_xgb <- auc_binary(yte, p_xgb)
      
      if (is.finite(auc_rf))  delta_rf  <- c(delta_rf,  auc_rf  - auc0)
      if (is.finite(auc_xgb)) delta_xgb <- c(delta_xgb, auc_xgb - auc0)
    }
  }
  
  c(
    rf  = mean(delta_rf,  na.rm = TRUE),
    xgb = mean(delta_xgb, na.rm = TRUE),
    n_rf  = sum(is.finite(delta_rf)),
    n_xgb = sum(is.finite(delta_xgb))
  )
}

perm_p_delta_auc <- function(y_pheno, X_taxa, age_vec, k, repeats, seed,
                             B = 200,  # increase to 1000 for final
                             within_strata = FALSE) {
  
  obs <- cv_mean_delta_auc(y_pheno, X_taxa, age_vec, k, repeats, seed)
  obs_rf  <- obs[["rf"]]
  obs_xgb <- obs[["xgb"]]
  
  null_rf  <- numeric(B)
  null_xgb <- numeric(B)
  
  set.seed(seed + 9999)
  for (b in seq_len(B)) {
    if (!within_strata) {
      y_perm <- sample(y_pheno)   # permute labels
    } else {
      # optional: permute within age strata (e.g., quartiles) if you want exchangeability conditional on age
      strata <- cut(age_vec, breaks = quantile(age_vec, probs = seq(0,1,0.25), na.rm = TRUE), include.lowest = TRUE)
      y_perm <- y_pheno
      for (s in levels(strata)) {
        ii <- which(strata == s)
        y_perm[ii] <- sample(y_perm[ii])
      }
    }
    
    out <- cv_mean_delta_auc(y_perm, X_taxa, age_vec, k, repeats, seed + 100*b)
    null_rf[b]  <- out[["rf"]]
    null_xgb[b] <- out[["xgb"]]
  }
  
  # one-sided: P(null >= observed)
  p_rf  <- (1 + sum(null_rf  >= obs_rf,  na.rm = TRUE)) / (1 + B)
  p_xgb <- (1 + sum(null_xgb >= obs_xgb, na.rm = TRUE)) / (1 + B)
  
  list(
    obs = obs,
    p_rf = p_rf,
    p_xgb = p_xgb,
    null_rf = null_rf,
    null_xgb = null_xgb
  )
}

run_one_bioproject_binary <- function(df_bp,
                                      k = 5,
                                      repeats = 20,
                                      seed = 1,
                                      min_n = 20,
                                      min_class = 10,
                                      do_age_task = TRUE) {
  
  df_bp <- df_bp %>% filter(!is.na(Phenotype), !is.na(Age))
  if (nrow(df_bp) < min_n) return(NULL)
  
  y_pheno <- factor(df_bp$Phenotype)
  if (nlevels(y_pheno) != 2) return(NULL)              # ONLY binary projects
  if (any(table(y_pheno) < min_class)) return(NULL)    # enough per class
  
  y_ageg <- make_age_group_median(df_bp$Age)
  if (do_age_task && any(table(y_ageg) < min_class)) do_age_task <- FALSE
  
  # Features (taxa only; age will be appended explicitly)
  meta_cols <- c("Age", "Phenotype", "BioProject")
  df_bp2 <- drop_allzero_features(df_bp, meta_cols = meta_cols)
  
  feat_cols <- setdiff(colnames(df_bp2), meta_cols)
  X_taxa <- as.matrix(df_bp2[, feat_cols, drop = FALSE])
  storage.mode(X_taxa) <- "numeric"
  
  age_vec <- as.numeric(df_bp2$Age)
  bp_id <- unique(df_bp2$BioProject)
  
  perm <- perm_p_delta_auc(y_pheno, X_taxa, age_vec, k, repeats, seed, B = 50)
  
  # Containers over ALL folds across repeats (Phenotype)
  auc_null_pheno     <- c()
  auc_full_rf_pheno  <- c()
  auc_full_xgb_pheno <- c()
  delta_rf_pheno     <- c()
  delta_xgb_pheno    <- c()
  
  # Importance accumulators over CV fits (Phenotype full models only)
  imp_rf_pheno_sum  <- NULL; imp_rf_pheno_n  <- 0L
  imp_xgb_pheno_sum <- NULL; imp_xgb_pheno_n <- 0L
  
  # Optional: AgeGroup task containers + importance
  auc_rf_ageg  <- c()
  auc_xgb_ageg <- c()
  imp_rf_age_sum  <- NULL; imp_rf_age_n  <- 0L
  imp_xgb_age_sum <- NULL; imp_xgb_age_n <- 0L
  
  for (r in seq_len(repeats)) {
    
    folds_pheno <- make_folds(y_pheno, k = k, seed = seed + 1000*r)
    if (do_age_task) folds_age <- make_folds(y_ageg, k = k, seed = seed + 2000*r)
    
    # ----- Phenotype folds -----
    for (i in seq_len(k)) {
      te <- folds_pheno[[i]]
      tr <- setdiff(seq_len(nrow(df_bp2)), te)
      
      ytr <- droplevels(y_pheno[tr])
      yte <- droplevels(y_pheno[te])
      
      # Guard against degenerate folds
      if (nlevels(ytr) != 2 || nlevels(yte) != 2) next
      if (any(table(ytr) < 2) || any(table(yte) < 1)) next
      
      Xtr_taxa <- X_taxa[tr, , drop = FALSE]
      Xte_taxa <- X_taxa[te, , drop = FALSE]
      age_tr   <- age_vec[tr]
      age_te   <- age_vec[te]
      
      # ---------------- Null: Age-only (logistic) ----------------
      # Fit with explicit level handling: positive = levels(ytr)[2]
      lev <- levels(ytr)
      ytr_bin <- ifelse(ytr == lev[2], 1, 0)
      
      fit_null <- glm(ytr_bin ~ age_tr, family = binomial())
      p_null <- as.numeric(predict(fit_null, newdata = data.frame(age_tr = age_te), type = "response"))
      auc0 <- auc_binary(yte, p_null)
      
      auc_null_pheno <- c(auc_null_pheno, auc0)
      
      # ---------------- Full RF: Taxa + Age ----------------
      Xtr_full <- cbind(Xtr_taxa, Age = age_tr)
      Xte_full <- cbind(Xte_taxa, Age = age_te)
      
      rf <- ranger::ranger(
        dependent.variable.name = "y",
        data = data.frame(y = ytr, Xtr_full, check.names = FALSE),
        probability = TRUE,
        num.trees = 1000,
        importance = "impurity",
        seed = seed + 10*r + i
      )
      
      prf <- as.matrix(predict(rf, data = data.frame(Xte_full, check.names = FALSE))$predictions)
      # probability for positive class = lev[2]
      p_rf <- if (lev[2] %in% colnames(prf)) prf[, lev[2]] else prf[, 2]
      auc_rf <- auc_binary(yte, p_rf)
      
      auc_full_rf_pheno <- c(auc_full_rf_pheno, auc_rf)
      delta_rf_pheno    <- c(delta_rf_pheno, auc_rf - auc0)
      
      imp_rf_pheno_sum <- add_named(imp_rf_pheno_sum, rf_imp_vec(rf))
      imp_rf_pheno_n <- imp_rf_pheno_n + 1L
      
      # ---------------- Full XGB: Taxa + Age ----------------
      dtr <- xgb.DMatrix(data = as.matrix(Xtr_full), label = ytr_bin)
      dte <- xgb.DMatrix(data = as.matrix(Xte_full))
      
      xgb_fit <- xgb.train(
        params = list(
          objective = "binary:logistic",
          eval_metric = "auc",
          max_depth = 6,
          eta = 0.05,
          subsample = 0.8,
          colsample_bytree = 0.8
        ),
        data = dtr,
        nrounds = 500,
        verbose = 0
      )
      
      p_xgb <- as.numeric(predict(xgb_fit, dte))
      auc_xgb <- auc_binary(yte, p_xgb)
      
      auc_full_xgb_pheno <- c(auc_full_xgb_pheno, auc_xgb)
      delta_xgb_pheno    <- c(delta_xgb_pheno, auc_xgb - auc0)
      
      imp_xgb_pheno_sum <- add_named(imp_xgb_pheno_sum, xgb_imp_vec(xgb_fit, colnames(Xtr_full)))
      imp_xgb_pheno_n <- imp_xgb_pheno_n + 1L
    }
    
    # ----- AgeGroup task (optional) -----
    if (do_age_task) {
      for (i in seq_len(k)) {
        te <- folds_age[[i]]
        tr <- setdiff(seq_len(nrow(df_bp2)), te)
        
        ytr <- droplevels(y_ageg[tr])
        yte <- droplevels(y_ageg[te])
        if (nlevels(ytr) != 2 || nlevels(yte) != 2) next
        
        Xtr_taxa <- X_taxa[tr, , drop = FALSE]
        Xte_taxa <- X_taxa[te, , drop = FALSE]
        
        # RF
        rf <- ranger::ranger(
          dependent.variable.name = "y",
          data = data.frame(y = ytr, Xtr_taxa, check.names = FALSE),
          probability = TRUE,
          num.trees = 1000,
          importance = "impurity",
          seed = seed + 30*r + i
        )
        prf <- as.matrix(predict(rf, data = data.frame(Xte_taxa, check.names = FALSE))$predictions)
        p_rf <- if ("Older" %in% colnames(prf)) prf[, "Older"] else prf[, 2]
        auc_rf_ageg <- c(auc_rf_ageg, auc_binary(yte, p_rf))
        imp_rf_age_sum <- add_named(imp_rf_age_sum, rf_imp_vec(rf))
        imp_rf_age_n <- imp_rf_age_n + 1L
        
        # XGB
        ytr_bin <- ifelse(ytr == "Older", 1, 0)
        dtr <- xgb.DMatrix(data = as.matrix(Xtr_taxa), label = ytr_bin)
        dte <- xgb.DMatrix(data = as.matrix(Xte_taxa))
        
        xgb_fit <- xgb.train(
          params = list(
            objective = "binary:logistic",
            eval_metric = "auc",
            max_depth = 6,
            eta = 0.05,
            subsample = 0.8,
            colsample_bytree = 0.8
          ),
          data = dtr,
          nrounds = 500,
          verbose = 0
        )
        p_xgb <- as.numeric(predict(xgb_fit, dte))
        auc_xgb_ageg <- c(auc_xgb_ageg, auc_binary(yte, p_xgb))
        imp_xgb_age_sum <- add_named(imp_xgb_age_sum, xgb_imp_vec(xgb_fit, colnames(Xtr_taxa)))
        imp_xgb_age_n <- imp_xgb_age_n + 1L
      }
    }
  }
  
  # Average importances across CV fits
  imp_rf_pheno <- if (imp_rf_pheno_n > 0) sort(imp_rf_pheno_sum / imp_rf_pheno_n, decreasing = TRUE) else numeric()
  imp_xgb_pheno <- if (imp_xgb_pheno_n > 0) sort(imp_xgb_pheno_sum / imp_xgb_pheno_n, decreasing = TRUE) else numeric()
  
  imp_rf_age <- if (imp_rf_age_n > 0) sort(imp_rf_age_sum / imp_rf_age_n, decreasing = TRUE) else numeric()
  imp_xgb_age <- if (imp_xgb_age_n > 0) sort(imp_xgb_age_sum / imp_xgb_age_n, decreasing = TRUE) else numeric()
  
  metrics <- data.frame(
    BioProject = bp_id,
    n = nrow(df_bp2),
    k = k,
    repeats = repeats,
    n_models_pheno = length(auc_null_pheno),
    
    pheno_auc_null_mean     = mean(auc_null_pheno,     na.rm = TRUE),
    pheno_auc_null_sd       = sd(auc_null_pheno,       na.rm = TRUE),
    
    pheno_auc_full_rf_mean  = mean(auc_full_rf_pheno,  na.rm = TRUE),
    pheno_auc_full_rf_sd    = sd(auc_full_rf_pheno,    na.rm = TRUE),
    
    pheno_auc_full_xgb_mean = mean(auc_full_xgb_pheno, na.rm = TRUE),
    pheno_auc_full_xgb_sd   = sd(auc_full_xgb_pheno,   na.rm = TRUE),
    
    pheno_delta_auc_rf_mean  = mean(delta_rf_pheno,    na.rm = TRUE),
    pheno_delta_auc_rf_sd    = sd(delta_rf_pheno,      na.rm = TRUE),
    
    pheno_delta_auc_xgb_mean = mean(delta_xgb_pheno,   na.rm = TRUE),
    pheno_delta_auc_xgb_sd   = sd(delta_xgb_pheno,     na.rm = TRUE),
    
    pheno_delta_auc_rf_perm_p  = perm$p_rf,
    pheno_delta_auc_xgb_perm_p = perm$p_xgb,
    
    do_age_task = do_age_task
  )
  
  if (do_age_task) {
    metrics$age_auc_rf_mean  <- mean(auc_rf_ageg,  na.rm = TRUE)
    metrics$age_auc_rf_sd    <- sd(auc_rf_ageg,    na.rm = TRUE)
    metrics$age_auc_xgb_mean <- mean(auc_xgb_ageg, na.rm = TRUE)
    metrics$age_auc_xgb_sd   <- sd(auc_xgb_ageg,   na.rm = TRUE)
    metrics$n_models_age     <- length(auc_rf_ageg)
  }
  
  list(
    metrics = metrics,
    importance = list(
      phenotype_rf  = data.frame(feature = names(imp_rf_pheno),  importance = as.numeric(imp_rf_pheno),  row.names = NULL),
      phenotype_xgb = data.frame(feature = names(imp_xgb_pheno), importance = as.numeric(imp_xgb_pheno), row.names = NULL),
      age_rf        = data.frame(feature = names(imp_rf_age),    importance = as.numeric(imp_rf_age),    row.names = NULL),
      age_xgb       = data.frame(feature = names(imp_xgb_age),   importance = as.numeric(imp_xgb_age),   row.names = NULL)
    )
  )
}

# ---------- main loop over datasets in a folder ----------

run_folder <- function(folder,
                       pattern = "\\.csv$",
                       k = 5,
                       repeats = 20,
                       seed = 1,
                       min_n = 20,
                       min_class = 10,
                       out_dir = file.path(folder, "model_outputs_binary")) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No CSV files found in folder: ", folder)
  
  all_metrics <- list()
  
  for (f in files) {
    message("Reading: ", f)
    dat <- read.csv(f, row.names = 1, header = TRUE, check.names = FALSE)
    
    req <- c("Age", "Phenotype", "BioProject")
    if (!all(req %in% colnames(dat))) {
      stop("Missing required columns in ", basename(f), ": ",
           paste(setdiff(req, colnames(dat)), collapse = ", "))
    }
    
    dat$BioProject <- as.character(dat$BioProject)
    dat$Phenotype  <- as.character(dat$Phenotype)
    dat$Age        <- as.numeric(dat$Age)
    
    for (bp in sort(unique(dat$BioProject))) {
      df_bp <- dat[dat$BioProject == bp, , drop = FALSE]
      
      res <- run_one_bioproject_binary(
        df_bp,
        k = k,
        repeats = repeats,
        seed = seed,
        min_n = min_n,
        min_class = min_class,
        do_age_task = TRUE
      )
      if (is.null(res)) next
      
      res$metrics$dataset <- basename(f)
      key <- paste(basename(f), bp, sep = " :: ")
      all_metrics[[key]] <- res$metrics
      
      prefix <- paste0(tools::file_path_sans_ext(basename(f)), "__", bp)
      
      write.csv(head(res$importance$phenotype_rf, 50),
                file = file.path(out_dir, paste0(prefix, "__imp_phenotype_full_rf_top50.csv")),
                row.names = FALSE)
      write.csv(head(res$importance$phenotype_xgb, 50),
                file = file.path(out_dir, paste0(prefix, "__imp_phenotype_full_xgb_top50.csv")),
                row.names = FALSE)
      
      if (isTRUE(res$metrics$do_age_task)) {
        write.csv(head(res$importance$age_rf, 50),
                  file = file.path(out_dir, paste0(prefix, "__imp_age_rf_top50.csv")),
                  row.names = FALSE)
        write.csv(head(res$importance$age_xgb, 50),
                  file = file.path(out_dir, paste0(prefix, "__imp_age_xgb_top50.csv")),
                  row.names = FALSE)
      }
      
      message("  Done BioProject: ", bp,
              " | Null(Age) AUC=", round(res$metrics$pheno_auc_null_mean, 3),
              " | Full RF AUC=", round(res$metrics$pheno_auc_full_rf_mean, 3),
              " Δ=", round(res$metrics$pheno_delta_auc_rf_mean, 3),
              " | Full XGB AUC=", round(res$metrics$pheno_auc_full_xgb_mean, 3),
              " Δ=", round(res$metrics$pheno_delta_auc_xgb_mean, 3))
    }
  }
  
  metrics_df <- dplyr::bind_rows(all_metrics)
  write.csv(metrics_df, file = file.path(out_dir, "all_bioproject_metrics_binary.csv"), row.names = FALSE)
  metrics_df
}

# ---------- Example run ----------
results <- run_folder(
  folder = "../data/clean",
  k = 5,
  repeats = 50,
  seed = 345,
  min_n = 20,
  min_class = 10
)

write.csv(results, file = "../results/auroc.csv", row.names = FALSE)

