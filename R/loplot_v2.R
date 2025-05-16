#' @title Loplot
#' @description Loplot, a visible figure of signal correction loplot
#' provide the visible figure of QC-based correction. See the details at the following References.
#' @param x the file before QC-RLS correction.
#' @param z the file after QC-RLS correction.
#' @param i a index for the name of variable.
#' @param fig_type either png or pdf. Defaults to pdf in shiftCor_v2.
#' @usage loplot(x,z,i,fig_type = 'pdf')
#' @export
#' @references statTarget: a streamlined tool for signal drift correction
#' and interpretations of quantitative mass spectrometry-based omics data.
#' Luan H, Ji F, Chen Y, Cai Z. 2018, Analytica Chimica Acta.
loplot_rlsc <- function(x, z, i, fig_type = 'pdf',
                        batch_wise = FALSE,
                        batch = NULL) {
  # - x, z: numeric matrices [features × injections]
  # - batch_wise: whether to do per-batch smoothing
  # - batch: a factor or character vector, length = ncol(x), giving each injection's batch
  
  # … all of your device + layout code stays the same …
  
  # pick up our batch info
  if (batch_wise) {
    if (is.null(batch) || length(batch) != ncol(x))
      stop("You must supply a `batch` vector of length ncol(x).")
    batches <- unique(batch)
    palette <- scales::hue_pal()(length(batches))
  }
  
  numY <- seq_len(ncol(x))
  qcid <- grep("QC", colnames(x))
  
  ##### RAW PANEL #####
  plot(numY, x[i,], pch=19, col="#F5C710",
       ylab="Intensity", xlab="Injection Order", main="Raw Peak")
  points(qcid, x[i, qcid], pch=19, col="blue")
  
  if (!batch_wise) {
    loe <- loess(x[i, qcid] ~ qcid)
    lines(numY, predict(loe, numY), col=rgb(0,0,0,0.3), lwd=4)
    legend("top", c("Sample","QC"), col=c("#F5C710","blue"),
           pch=19, bty="n", cex=0.75, horiz=TRUE)
  } else {
    legend_labels <- c("QC")
    legend_cols   <- "blue"
    for (j in seq_along(batches)) {
      idx   <- which(batch == batches[j])
      qc_b  <- intersect(qcid, idx)
      if (length(qc_b) < 5) next
      loe_b <- loess(x[i, qc_b] ~ qc_b)
      lines(qc_b, predict(loe_b, qc_b), col=palette[j], lwd=2)
      legend_labels <- c(legend_labels, paste0("Batch ", batches[j]))
      legend_cols   <- c(legend_cols, palette[j])
    }
    legend("topright", legend_labels, col=legend_cols,
           pch=c(19, rep(NA, length(batches))), lty=c(NA, rep(1, length(batches))),
           bty="n", cex=0.7)
  }
  
  ##### CORRECTED PANEL #####
  plot(numY, z[i,], pch=19, col="#F5C710",
       ylab="Intensity", xlab="Injection Order", main="Corrected Peak")
  points(qcid, z[i, qcid], pch=19, col="blue")
  
  if (!batch_wise) {
    loe_n <- loess(z[i, qcid] ~ qcid)
    lines(numY, predict(loe_n, numY), col=rgb(0,0,0,0.3), lwd=4)
    legend("top", c("Sample","QC"), col=c("#F5C710","blue"),
           pch=19, bty="n", cex=0.75, horiz=TRUE)
  } else {
    legend_labels <- c("QC")
    legend_cols   <- "blue"
    for (j in seq_along(batches)) {
      idx   <- which(batch == batches[j])
      qc_b  <- intersect(qcid, idx)
      if (length(qc_b) < 5) next
      loe_bc <- loess(z[i, qc_b] ~ qc_b)
      lines(qc_b, predict(loe_bc, qc_b), col=palette[j], lwd=2, lty=2)
      legend_labels <- c(legend_labels, paste0("Batch ", batches[j]))
      legend_cols   <- c(legend_cols, palette[j])
    }
    legend("topright", legend_labels, col=legend_cols,
           pch=c(19, rep(NA, length(batches))), lty=c(NA, rep(2, length(batches))),
           bty="n", cex=0.7)
  }
  
  dev.off()
}