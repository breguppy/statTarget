#' @name QCRFSC_oob
#' @title QC‑RFSC with OOB extraction
#' @description
#'   Performs QC‑based RF signal correction and returns both the corrected
#'   data and the per‑metabolite out‑of‑bag MSEs.
#' @param samPeno   Path to the sample‑list CSV (with columns: sample, class, batch, order, …).
#' @param samFile   Path to the profile CSV (first row = header of sample names).
#' @param Frule     Fractional rule for filtering features (default 0.8).
#' @param ntree     Number of trees for each RF (default 500).
#' @param seed      Optional random seed (default NULL).
#' @param imputeM   Imputation method: "KNN", "min", "minHalf", or "median" (default "KNN").
#' @return A list with two elements:
#'   * corrected_matrix: a matrix [metabolites × samples] of drift‑corrected intensities  
#'   * oob_mse_per_metabolite: numeric vector of length = #metabolites
#' @export
QCRFSC_oob <- function(
    samPeno,
    samFile,
    Frule   = 0.8,
    ntree   = 500,
    seed    = NULL,
    imputeM = "KNN"
) {
  #— load deps
  library(plyr)            # for arrange()
  library(impute)          # for impute.knn()
  library(randomForest)    # for randomForest()
  
  #— 1) Read & align metadata + profile data
  samPeno_df <- read.csv(samPeno, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  samPeno_df <- arrange(samPeno_df, order)
  
  prof_raw <- read.csv(samFile, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  prof_mat <- t(prof_raw)
  colnames(prof_mat) <- prof_mat[1, ]
  prof_df <- as.data.frame(prof_mat[-1, ], stringsAsFactors = FALSE)
  rownames(prof_df) <- prof_df$name
  
  # ensure samples line up
  idx <- match(samPeno_df$sample, prof_df$name)
  if (any(is.na(idx))) stop("Some samples in samPeno are missing from samFile")
  samFP <- prof_df[idx, ]
  
  #— 2) Build numeric matrix & replace zeros with NA
  mat <- as.matrix(samFP[, setdiff(colnames(samFP), "name")])
  mat[mat < 0] <- 0
  mat[mat == 0] <- NA
  
  #— 3) Filter features by Frule
  classF <- addNA(as.factor(samPeno_df$class))
  keep_cols <- sapply(seq_len(ncol(mat)), function(j) {
    freqs <- tapply(mat[, j], classF, function(x) sum(is.na(x)) / length(x))
    any(freqs <= (1 - Frule))
  })
  mat <- mat[, keep_cols, drop = FALSE]
  
  #— 4) Impute
  if (imputeM == "KNN") {
    inputedData <- impute.knn(mat, rowmax = 0.99, colmax = 0.99, maxp = 15000)$data
  } else if (imputeM == "min") {
    minValue <- function(x, grp) {
      g <- as.factor(as.numeric(grp))
      for (i in seq_len(nrow(x)))
        for (j in seq_len(ncol(x)))
          if (is.na(x[i,j]))
            x[i,j] <- tapply(x[,j], g, min, na.rm = TRUE)[g[i]]
      x
    }
    inputedData <- minValue(mat, classF)
  } else if (imputeM == "minHalf") {
    minHalfValue <- function(x, grp) {
      g <- as.factor(as.numeric(grp))
      for (i in seq_len(nrow(x)))
        for (j in seq_len(ncol(x)))
          if (is.na(x[i,j]))
            x[i,j] <- tapply(x[,j], g, min, na.rm = TRUE)[g[i]] / 2
      x
    }
    inputedData <- minHalfValue(mat, classF)
  } else if (imputeM == "median") {
    medianValue <- function(x, grp) {
      g <- as.factor(as.numeric(grp))
      for (i in seq_len(nrow(x)))
        for (j in seq_len(ncol(x)))
          if (is.na(x[i,j]))
            x[i,j] <- tapply(x[,j], g, median, na.rm = TRUE)[g[i]]
      x
    }
    inputedData <- medianValue(mat, classF)
  } else {
    stop("Unknown imputation method: ", imputeM)
  }
  
  #— 5) Prepare for QC‑RFSC
  #    rows = metabolites, cols = samples
  dat  <- t(inputedData)
  order_positions <- seq_len(ncol(dat))
  
  #— 6) Fit RF per metabolite, collect OOB, correct drift
  REGfit <- function(x, y, ntree, seed) {
    n_met    <- nrow(x)
    oob_errs <- numeric(n_met)
    if (!is.null(seed)) set.seed(seed)
    qc_idx <- grep("QC", colnames(x))
    
    for (i in seq_len(n_met)) {
      rf_mod <- randomForest(
        x     = data.frame(order = qc_idx),
        y     = as.numeric(x[i, qc_idx]),
        ntree = ntree
      )
      # store final OOB MSE
      oob_errs[i] <- rf_mod$mse[length(rf_mod$mse)]
      # predict intensities at all positions
      preds      <- predict(rf_mod, newdata = data.frame(order = y))
      x[i, ]     <- x[i, ] / preds
    }
    list(corrected = x, oob_mse = oob_errs)
  }
  
  fit_res  <- REGfit(dat, order_positions, ntree = ntree, seed = seed)
  loessDat <- fit_res$corrected
  oob_vec  <- fit_res$oob_mse
  
  #— 7) Restore dimnames
  if (!is.null(rownames(dat))) rownames(loessDat) <- rownames(dat)
  if (!is.null(colnames(dat))) colnames(loessDat) <- colnames(dat)
  
  #— 8) Return
  list(
    corrected_matrix       = loessDat,
    oob_mse_per_metabolite = oob_vec
  )
}

