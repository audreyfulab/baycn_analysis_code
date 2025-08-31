########### read performance evaluation ###########
# G2
g2_all <- read.delim("g2_prerec_mse2_all_methods_table.tsv", header = TRUE, sep = "\t", na.strings = "-")

# pc
pc_all <- read.delim("pc_prerec_mse2_all_methods_table.tsv", header = TRUE, sep = "\t", na.strings = "-")

# nc11
nc11_all <- read.delim("nc11_prerec_mse2_all_methods_table_merged.tsv", header = TRUE, sep = "\t", na.strings = "-")


########### visualization ###########
methods <- c(setdiff(unique(nc11_all$Method), "baycn"), "baycn")
# "baycn", "mc3", "order", "partition", "bcdag", "scanBMA"
mycolors <- c("orange", "green", "magenta", "blue", "brown", "black")
mycolors.lines <- adjustcolor(mycolors, alpha.f = 0.5)

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

plot_with_SD_bars <- function(data, xvar, yvar, xsdvar, ysdvar, 
                                 method_col = "Method", methods = NULL, 
                                 colors = NULL, pch = 16, cex = 1.2, 
                                 xlab = NULL, ylab = NULL, 
                                 legend=FALSE, legend_position = NULL,
                                 diagonal=FALSE, 
                                 sd_colors = NULL, sd_lty = 1, sd_lwd=0.8,...) {
  
  # Map each method to color
  method_colors <- setNames(colors, methods)
  line_colors <- setNames(sd_colors, methods)
  
  # Plot base
  plot(data[[xvar]], data[[yvar]], type = "n",
       xlab = xlab %||% xvar, ylab = ylab %||% yvar, ...)
  if (diagonal){
    abline(a = 0, b = 1, lty=3)
  }
  
  for (m in methods) {
    idx <- which(data[[method_col]] == m)
    x <- data[[xvar]][idx]
    y <- data[[yvar]][idx]
    xsd <- data[[xsdvar]][idx]
    ysd <- data[[ysdvar]][idx]
    col <- method_colors[[m]]
    sd_col <- line_colors[[m]]
    
    # Error bars
    arrows(x - xsd, y, x + xsd, y, angle = 90, code = 3, length = 0.02, col = sd_col, lty = sd_lty, lwd = sd_lwd)
    arrows(x, y - ysd, x, y + ysd, angle = 90, code = 3, length = 0.02, col = sd_col, lty = sd_lty, lwd = sd_lwd)
    
    # Points
    points(x, y, col = col, pch = pch, cex = cex)
  }
  
  if (legend){
    legend(legend_position,
           legend = methods,
           col = colors,
           pch = pch,
           pt.cex = cex,
           bty = "n",  # no box around legend
           cex = 0.9,
           y.intersp = 0.8)
  }
  
}


# Adjust spacing
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

pdf("Plot_power_precision_g2_nc11_pc.pdf")
par(mar = c(3.5, 3.5, 0.5, 0.5), mgp = c(2, 0.5, 0))
par(mfrow=c(3,3))

#### g2
# recall vs precision
plot_with_SD_bars(g2_all, 
                  xvar = "Precision.Mean", yvar = "Recall.Mean", 
                  xsdvar = "Precision.SD", ysdvar = "Recall.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  legend = TRUE, legend_position = "topleft",
                  xlab = "Precision", ylab = "Power", diagonal = TRUE, xlim=c(0, 1), ylim=c(0,1))

# MSE2 vs precision
plot_with_SD_bars(g2_all, 
                  xvar = "Precision.Mean", yvar = "MSE.Mean", 
                  xsdvar = "Precision.SD", ysdvar = "MSE.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  xlab = "Precision", ylab = "MSE2", xlim=c(0, 1), ylim=c(0, 0.3))

# MSE2 vs recall
plot_with_SD_bars(g2_all, 
                  xvar = "Recall.Mean", yvar = "MSE.Mean", 
                  xsdvar = "Recall.SD", ysdvar = "MSE.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  xlab = "Power", ylab = "MSE2", xlim=c(0, 1), ylim=c(0, 0.3))

#### NC11
# recall vs precision
plot_with_SD_bars(nc11_all, 
                     xvar = "Precision.Mean", yvar = "Recall.Mean", 
                     xsdvar = "Precision.SD", ysdvar = "Recall.SD",
                     methods = methods,
                     colors = mycolors,
                     sd_colors = mycolors.lines,
                     sd_lty = 1,
                     xlab = "Precision", ylab = "Power", 
                     diagonal = TRUE, xlim=c(0, 1), ylim=c(0,1))

# MSE2 vs precision
plot_with_SD_bars(nc11_all, 
                  xvar = "Precision.Mean", yvar = "MSE.Mean", 
                  xsdvar = "Precision.SD", ysdvar = "MSE.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  sd_lty = 1,
                  xlab = "Precision", ylab = "MSE2", 
                  legend_position = 'topleft',
                  diagonal = FALSE, xlim=c(0, 1), ylim=c(0,0.1))

# MSE2 vs recall
plot_with_SD_bars(nc11_all, 
                  xvar = "Recall.Mean", yvar = "MSE.Mean", 
                  xsdvar = "Recall.SD", ysdvar = "MSE.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  sd_lty = 1,
                  xlab = "Power", ylab = "MSE2", 
                  legend_position = 'topleft',
                  xlim=c(0, 1), ylim=c(0, 0.1))


#### pc
# recall vs precision
plot_with_SD_bars(pc_all, 
                  xvar = "Precision.Mean", yvar = "Recall.Mean", 
                  xsdvar = "Precision.SD", ysdvar = "Recall.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  xlab = "Precision", ylab = "Power", diagonal = TRUE, xlim=c(0, 1), ylim=c(0,1))

# MSE2 vs precision
plot_with_SD_bars(pc_all, 
                  xvar = "Precision.Mean", yvar = "MSE.Mean", 
                  xsdvar = "Precision.SD", ysdvar = "MSE.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  xlab = "Precision", ylab = "MSE2", xlim=c(0, 1), ylim=c(0, 0.14))

# MSE2 vs recall
plot_with_SD_bars(pc_all, 
                  xvar = "Recall.Mean", yvar = "MSE.Mean", 
                  xsdvar = "Recall.SD", ysdvar = "MSE.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  xlab = "Power", ylab = "MSE2", xlim=c(0, 1), ylim=c(0, 0.14))

dev.off()

pdf ("Plot_power_precision_g2_fully_connected.pdf")
plot_with_SD_bars(g2_all, 
                  xvar = "Precision.Mean", yvar = "Recall.Mean", 
                  xsdvar = "Precision.SD", ysdvar = "Recall.SD",
                  methods = methods,
                  colors = mycolors,
                  sd_colors = mycolors.lines,
                  legend = TRUE, legend_position = "topleft",
                  xlab = "Precision", ylab = "Power", diagonal = TRUE, xlim=c(0, 1), ylim=c(0,1))
dev.off()
