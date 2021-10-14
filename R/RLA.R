rlaPlots<-function (inputdata, type = c("ag", "wg"), cols = NULL, cex.axis = 0.8, 
                    las = 2, ylim = c(-2, 2), oma = c(7, 4, 4, 2) + 0.1, ...) 
{
  type <- match.arg(type)
  groups <- factor(inputdata[, 1], levels = unique(inputdata[,1]))
  unique.groups <- levels(groups)
  
  if (is.null(cols)) 
  cols <- rainbow(length(unique.groups))
  box_cols <- c(rep(NA, length(rownames(inputdata))))
  inputdata<-inputdata[order(inputdata[,1]),]
  for (ii in 1:length(inputdata[, 1])) {
    
    box_cols[ii] <- cols[match(inputdata[, 1][ii],unique.groups)]
  }
  
  if (type == "wg") {
    out_data <- data.frame()
    for (grp in unique.groups) {
      submat <- inputdata[which(inputdata[, 1] == grp), 
                          -1]
      med_vals <- apply(submat, 2, median)
      swept_mat <- sweep(submat, 2, med_vals, "-")
      out_data <- rbind(out_data, swept_mat)
    }
  }
  else {
    med_vals <- apply(inputdata[, -1], 2, median)
    out_data <- sweep(inputdata[, -1], 2, med_vals, "-")
  }
  # ylim<-round(range(out_data))*0.9 # control range. 
  boxplot(t(out_data), cex.axis = cex.axis, las = las, col=box_cols,
          oma = oma, outline = FALSE, ...)
  abline(h = 0)
}


