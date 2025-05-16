#' @title Loplot (batch-aware)
#' @description Visible QC-based LOESS correction plots, with an option for batch-wise smoothing.
#' @param x          matrix of raw intensities (features × injections)
#' @param z          matrix of corrected intensities (features × injections)
#' @param i          row index of the feature to plot
#' @param fig_type   either "png" or "pdf"
#' @param batch_wise logical; if TRUE, draw one grey LOESS curve per batch
#' @param batch      vector (length = ncol(x)) of batch labels for each injection
#' @export
loplot_rlsc <- function(x, z, i,
                        fig_type   = 'pdf',
                        batch_wise = FALSE,
                        batch      = NULL) {
  # sanity check
  if (batch_wise) {
    if (is.null(batch) || length(batch) != ncol(x))
      stop("When batch_wise=TRUE you must supply a ‘batch’ vector of length ncol(x).")
    batches <- unique(batch)
  }
  
  # set up output folder
  met_name <- gsub("[^[:alnum:]]", "_", rownames(x)[i])
  out_dir  <- file.path(getwd(),
                        "statTarget/shiftCor/After_shiftCor/loplot")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  
  # open device
  filename <- file.path(out_dir,
                        paste0(met_name, "_", i,
                               if (fig_type=="png") ".png" else ".pdf"))
  if (fig_type == 'png') {
    png(filename, width=2100, height=2100, res=300)
    par(mar = c(5,5,4,2) + 0.1)
    layout(matrix(1:2, nrow=2))
  } else {
    pdf(filename, width=6, height=6)
    layout(matrix(1:2, nrow=2))
  }
  
  numY <- seq_len(ncol(x))
  qcid <- grep("QC", colnames(x))
  
  #### RAW PANEL ####
  plot(numY, x[i,], pch=19, col="#F5C710",
       ylab="Intensity", xlab="Injection Order", main="Raw Peak")
  points(qcid, x[i, qcid], pch=19, col="blue")
  
  if (!batch_wise) {
    # single global grey LOESS
    loe <- loess(x[i, qcid] ~ qcid)
    lines(numY, predict(loe, numY),
          col=rgb(0,0,0,0.3), lwd=4)
  } else {
    # one grey LOESS curve per batch
    for (b in batches) {
      idx  <- which(batch == b)
      qc_b <- intersect(qcid, idx)
      if (length(qc_b) < 3) next
      loe_b <- loess(x[i, qc_b] ~ qc_b)
      lines(qc_b, predict(loe_b, qc_b),
            col=rgb(0,0,0,0.3), lwd=4)
    }
  }
  
  legend("top", c("Sample","QC"), col=c("#F5C710","blue"),
         pch=19, bty="n", cex=0.75, horiz=TRUE)
  
  #### CORRECTED PANEL ####
  plot(numY, z[i,], pch=19, col="#F5C710",
       ylab="Intensity", xlab="Injection Order", main="Corrected Peak")
  points(qcid, z[i, qcid], pch=19, col="blue")
  
  if (!batch_wise) {
    loe_n <- loess(z[i, qcid] ~ qcid)
    lines(numY, predict(loe_n, numY),
          col=rgb(0,0,0,0.3), lwd=4)
  } else {
    for (b in batches) {
      idx   <- which(batch == b)
      qc_b  <- intersect(qcid, idx)
      if (length(qc_b) < 3) next
      loe_bc <- loess(z[i, qc_b] ~ qc_b)
      lines(qc_b, predict(loe_bc, qc_b),
            col=rgb(0,0,0,0.3), lwd=4)
    }
  }
  
  legend("top", c("Sample","QC"), col=c("#F5C710","blue"),
         pch=19, bty="n", cex=0.75, horiz=TRUE)
  
  dev.off()
}
