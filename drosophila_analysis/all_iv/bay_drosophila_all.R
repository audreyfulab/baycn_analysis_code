cat("Script started at:", format(Sys.time()), "\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript bay_drosophila_all.R <iterations> <run_number> <thinTo>")
}
iterations <- as.numeric(args[1])
run_number <- args[2]
thinTo <- as.numeric(args[3])

cat("Iterations:", iterations, "\n")
cat("Run number:", run_number, "\n")
cat("ThinTo:", thinTo, "\n")

# Define home directory
home_dir <- "/wsu/home/ht/ht26/ht2699"

# Load libraries
library ("baycn", lib=file.path(home_dir, "Rpackages"))

# Use data from the package
dros_cont_reordered <- drosophila$continuous[,c(7:21, 1:6)]

# Adjacency matrices 
am_all <- matrix(1, nrow = 21, ncol = 21)
diag(am_all) <- 0
am_all[16:21,] <- 0
rownames(am_all) <- colnames(dros_cont_reordered)
colnames(am_all) <- rownames(am_all)

# Run baycn
# Three independent long runs
baycn_dros <- mhEdge(data = dros_cont_reordered,
                     adjMatrix = am_all,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     nCPh = 6,
                     nGV = 0,
                     pmr = FALSE,
                     burnIn = 0.2,
                     iterations = iterations,
                     thinTo = thinTo,
                     progress = TRUE)

save(baycn_dros, file=paste0("dros_baycn_", iterations, "_", run_number, ".RData"))

cat("Script completed at:", format(Sys.time()), "\n")
