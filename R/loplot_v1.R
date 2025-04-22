### Loplot, a visible figure of signal correction loplot provide the visible figure of QC-based
### correction.  x The file before QC-based correction.  z The file after QC-based correction.  i
### The name of variable.  The plot of RSD distribution and data table.
loplot_rfsc <- function(x, z, i) {
    # x is the Corrected Data
    cn <- colnames(x)
    qcid <- grep("QC", cn)
    meantem <- mean(z[i, qcid])
    sdtem <- sd(z[i, qcid])
    rsdcutoffU <- 0.3
    rsdcutoffL <- 0.15
    # rsdcutoffX <- 1.0
    sdLineU <- rsdcutoffU * meantem
    sdLineL <- rsdcutoffL * meantem
    cuttemU <- meantem + sdLineU
    cuttemD <- meantem - sdLineU
    cuttemUL <- meantem + sdLineL
    cuttemDL <- meantem - sdLineL
    
    meantemR <- mean(x[i, qcid])
    sdtemR <- sd(x[i, qcid])
    rsdcutoffRU <- 0.3
    rsdcutoffRL <- 0.15
    # rsdcutoffRX <- 1.0
    sdLineRU <- rsdcutoffRU * meantemR
    sdLineRL <- rsdcutoffRL * meantemR
    cuttemRU <- meantemR + sdLineRU
    cuttemRD <- meantemR - sdLineRU
    cuttemRUL <- meantemR + sdLineRL
    cuttemRDL <- meantemR - sdLineRL
    
    met_name <- gsub("[^[:alnum:]]", "_", rownames(x)[i])
    RSD30_CV=paste(met_name,"_", i,".png", sep="")
    #RSD30_CV = paste(rownames(x)[i], "_", i, ".pdf", sep = "")
    dirout.loplot <- paste(getwd(), "/statTarget/shiftCor/After_shiftCor/loplot", sep = "")
    if (!file.exists(dirout.loplot)){
      dir.create(dirout.loplot,showWarnings = FALSE)
    }
    png(
      filename = paste(dirout.loplot, RSD30_CV, sep = "/"),
      width = 2100,
      height = 2100,      
      res = 300
    )
    graphics::layout(matrix(1:2, nrow = 2))
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    numY <- 1:dim(x)[2]
    graphics::plot(numY, x[i, ], pch = 19, col = "#F5C710", ylab = c("Intensity"), xlab = c("Injection Order"), 
        main = "Raw Peak")
    points(qcid, x[i, qcid], pch = 19, col = "blue")
    abline(h = cuttemRU, col = "#D55E00", lwd = 1.2, lty = 1)
    abline(h = cuttemRD, col = "#D55E00", lwd = 1.2, lty = 1)
    abline(h = cuttemRUL, col = rgb(0, 0, 0, 0.3), lwd = 1.2, lty = 2)
    abline(h = cuttemRDL, col = rgb(0, 0, 0, 0.3), lwd = 1.2, lty = 2)
    mtext("+SD@30%CV", side = 4, line = 0.3, at = cuttemRU, col = "#D55E00", cex = 0.3, las = 1)
    mtext("-SD@30%CV", side = 4, line = 0.3, at = cuttemRD, col = "#D55E00", cex = 0.3, las = 1)
    mtext("+SD@15%CV", side = 4, line = 0.3, at = cuttemRUL, col = "#D55E00", cex = 0.3, las = 1)
    mtext("-SD@15%CV", side = 4, line = 0.3, at = cuttemRDL, col = "#D55E00", cex = 0.3, las = 1)
    
    legend("top", c("Sample", "QC"), col = c("#F5C710", "blue"), lty = 1, pch = 19, bty = "n", cex = 0.75, 
        horiz = TRUE)
    # lines(qcid,x[i,qcid],col=rgb(0,0,0,0.3),lwd=4) loe <- loess(x[i,qcid]~qcid)
    # points(numY,predict(loe,numY),type='l',col=rgb(0,0,0,0.3),lwd=4)
    
    graphics::plot(numY, z[i, ], pch = 19, col = "#F5C710", ylab = c("Intensity"), xlab = c("Injection Order"), 
        main = "Corrected Peak")
    points(qcid, z[i, qcid], pch = 19, col = "blue")
    abline(h = cuttemUL, col = rgb(0, 0, 0, 0.3), lwd = 1.2, lty = 2)
    abline(h = cuttemDL, col = rgb(0, 0, 0, 0.3), lwd = 1.2, lty = 2)
    abline(h = cuttemU, col = "#D55E00", lwd = 1.2, lty = 1)
    abline(h = cuttemD, col = "#D55E00", lwd = 1.2, lty = 1)
    mtext("+SD@30%CV", side = 4, line = 0.3, at = cuttemU, col = "#D55E00", cex = 0.3, las = 1)
    mtext("-SD@30%CV", side = 4, line = 0.3, at = cuttemD, col = "#D55E00", cex = 0.3, las = 1)
    mtext("+SD@15%CV", side = 4, line = 0.3, at = cuttemUL, col = "#D55E00", cex = 0.3, las = 1)
    mtext("-SD@15%CV", side = 4, line = 0.3, at = cuttemDL, col = "#D55E00", cex = 0.3, las = 1)
    
    # lines(qcid,z[i,qcid],col=rgb(0,0,0,0.3),lwd=4) loe_n <- loess(z[i,qcid]~qcid)
    # points(numY,predict(loe_n,numY),type='l',col=rgb(0,0,0,0.3),lwd=4)
    legend("top", c("Sample", "QC"), col = c("#F5C710", "blue"), lty = 1, pch = 19, bty = "n", cex = 0.75, 
        horiz = TRUE)
    dev.off()
}
