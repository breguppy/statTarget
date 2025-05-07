#' @name QCRFSC_cv_rmse
#' @title QC‑RFSC cross‑validated RMSE per metabolite
#' @description
#'   Reads sample metadata and profile data, applies QC‑RFSC drift correction,
#'   but instead of OOB error returns CV‑RMSE (leave‑one‑out or K‑fold) per metabolite.
#'
#' @param samPeno   Path to the sample‑list CSV (columns: sample, class, batch, order, ...).
#' @param samFile   Path to the profile CSV (first row = sample names header).
#' @param Frule     Fractional rule for filtering features (default 0.8).
#' @param ntree     Number of trees for each RF (default 500).
#' @param seed      Optional random seed (default NULL).
#' @param imputeM   Imputation: "KNN", "min", "minHalf", or "median" (default "KNN").
#' @param K         Number of folds for CV; if NULL or >= n_QC, uses leave‑one‑out (default NULL).
#'
#' @return Numeric vector of length = #metabolites: the CV‑RMSE for each metabolite.
#' @export
QCRFSC_cv_rmse <- function(
    samPeno,
    samFile,
    Frule   = 0.8,
    ntree   = 500,
    seed    = NULL,
    imputeM = "KNN",
    K       = NULL
) {
  #--- Dependencies
  library(plyr)         # for arrange()
  library(impute)       # for impute.knn()
  library(randomForest) # for randomForest()
  
  #--- 1) Read metadata & profile data
  samp <- read.csv(samPeno, header=TRUE, check.names=FALSE,
                   stringsAsFactors=FALSE)
  samp <- arrange(samp, order)
  
  prof_raw <- read.csv(samFile, header=FALSE,
                       check.names=FALSE, stringsAsFactors=FALSE)
  prof_mat <- t(prof_raw)
  colnames(prof_mat) <- prof_mat[1,]
  prof_df <- as.data.frame(prof_mat[-1,], stringsAsFactors=FALSE)
  rownames(prof_df) <- prof_df$name
  
  # align sample order
  idx <- match(samp$sample, prof_df$name)
  if(any(is.na(idx))) stop("Samples in metadata missing from profile file.")
  data_mat <- as.matrix(prof_df[idx, setdiff(colnames(prof_df), 'name')])
  
  #--- 2) Clean & replace zeros
  data_mat[data_mat < 0] <- 0
  data_mat[data_mat == 0] <- NA
  
  #--- 3) Filter features by Frule
  classF <- addNA(as.factor(samp$class))
  keep <- sapply(seq_len(ncol(data_mat)), function(j) {
    frq <- tapply(data_mat[,j], classF, function(x) sum(is.na(x))/length(x))
    any(frq <= (1 - Frule))
  })
  data_mat <- data_mat[, keep, drop=FALSE]
  
  #--- 4) Impute missing values
  if(imputeM == "KNN") {
    imp <- impute.knn(data_mat, rowmax=0.99, colmax=0.99, maxp=15000)
    imp_data <- imp$data
  } else if(imputeM == "min") {
    minVal <- function(x, g) {
      gr <- as.factor(as.numeric(g))
      for(i in seq_len(nrow(x))) for(j in seq_len(ncol(x))) {
        if(is.na(x[i,j])) x[i,j] <- tapply(x[,j], gr, min, na.rm=TRUE)[gr[i]]
      }
      x
    }
    imp_data <- minVal(data_mat, classF)
  } else if(imputeM == "minHalf") {
    minH <- function(x, g) {
      gr <- as.factor(as.numeric(g))
      for(i in seq_len(nrow(x))) for(j in seq_len(ncol(x))) {
        if(is.na(x[i,j])) x[i,j] <- tapply(x[,j], gr, min, na.rm=TRUE)[gr[i]]/2
      }
      x
    }
    imp_data <- minH(data_mat, classF)
  } else if(imputeM == "median") {
    medVal <- function(x, g) {
      gr <- as.factor(as.numeric(g))
      for(i in seq_len(nrow(x))) for(j in seq_len(ncol(x))) {
        if(is.na(x[i,j])) x[i,j] <- tapply(x[,j], gr, median, na.rm=TRUE)[gr[i]]
      }
      x
    }
    imp_data <- medVal(data_mat, classF)
  } else {
    stop("Unknown imputeM: ", imputeM)
  }
  
  #--- 5) Prepare for CV
  mat_t <- t(imp_data)  # metabolites × samples
  qc_idx <- grep("QC", colnames(mat_t))
  n_met <- nrow(mat_t)
  n_qc  <- length(qc_idx)
  
  #--- 6) Cross‑validated RMSE per metabolite
  if(!is.null(seed)) set.seed(seed)
  if(is.null(K) || K >= n_qc) {
    folds <- seq_len(n_qc)  # LOO
  } else {
    folds <- sample(rep(seq_len(K), length.out=n_qc))
  }
  
  rmse <- numeric(n_met)
  for(i in seq_len(n_met)) {
    y_qc <- as.numeric(mat_t[i, qc_idx])
    preds <- numeric(n_qc)
    for(f in unique(folds)) {
      test_i  <- which(folds == f)
      train_i <- setdiff(seq_len(n_qc), test_i)
      rf_mod <- randomForest(
        x     = data.frame(order=qc_idx[train_i]),
        y     = y_qc[train_i],
        ntree = ntree
      )
      preds[test_i] <- predict(
        rf_mod,
        newdata = data.frame(order=qc_idx[test_i])
      )
    }
    rmse[i] <- sqrt(mean((y_qc - preds)^2, na.rm=TRUE))
  }
  
  #--- 7) Return CV‑RMSE vector
  rmse
}
