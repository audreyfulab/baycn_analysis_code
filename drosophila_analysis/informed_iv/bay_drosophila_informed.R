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

# Adjacency matrices 
am_informed <- matrix(0, nrow = 21, ncol = 21)
rownames(am_informed) <- colnames(drosophila$continuous)
colnames(am_informed) <- rownames(am_informed)

am_informed['SM', 'Mef2_6_8h'] <- 1
am_informed['SM', 'Mef2_8_10h'] <- 1
am_informed['SM', 'Mef2_10_12h'] <- 1

am_informed['Meso_SM', 'Mef2_10_12h'] <- 1
am_informed['Meso_SM', 'Mef2_2_4h'] <- 1
am_informed['Meso_SM', 'Mef2_4_6h'] <- 1
am_informed['Meso_SM', 'Mef2_6_8h'] <- 1
am_informed['Meso_SM', 'Mef2_8_10h'] <- 1

am_informed['Meso_SM', 'Tin_2_4h'] <- 1
am_informed['Meso_SM', 'Tin_4_6h'] <- 1
am_informed['Meso_SM', 'Tin_6_8h'] <- 1

am_informed['VM_SM', 'Twi_2_4h'] <- 1
am_informed['VM_SM', 'Bin_8_10h'] <- 1
am_informed['VM_SM', 'Mef2_8_10h'] <- 1
am_informed['VM_SM', 'Bin_10_12h'] <- 1
am_informed['VM_SM', 'Mef2_10_12h'] <- 1

am_informed['Meso', 'Twi_2_4h'] <- 1

am_informed['VM', 'Bap_6_8h'] <- 1
am_informed['VM', 'Bin_6_8h'] <- 1
am_informed['VM', 'Bin_8_10h'] <- 1
am_informed['VM', 'Bin_10_12h'] <- 1

am_informed['CM', 'Bin_10_12h'] <- 1

am_informed['Mef2_2_4h', 'Mef2_10_12h'] <- 1
am_informed['Mef2_8_10h', 'Mef2_10_12h'] <- 1
am_informed['Mef2_2_4h', 'Mef2_6_8h'] <- 1
am_informed['Mef2_4_6h', 'Mef2_6_8h'] <- 1
am_informed['Mef2_6_8h', 'Mef2_8_10h'] <- 1

am_informed['Tin_2_4h', 'Tin_4_6h'] <- 1
am_informed['Tin_2_4h', 'Twi_4_6h'] <- 1
am_informed['Tin_4_6h', 'Twi_4_6h'] <- 1
am_informed['Tin_4_6h', 'Twi_6_8h'] <- 1
am_informed['Twi_2_4h', 'Twi_4_6h'] <- 1
am_informed['Twi_4_6h', 'Twi_6_8h'] <- 1
am_informed['Tin_4_6h', 'Mef2_2_4h'] <- 1
am_informed['Tin_4_6h', 'Mef2_4_6h'] <- 1
am_informed['Twi_4_6h', 'Mef2_4_6h'] <- 1
am_informed['Twi_6_8h', 'Mef2_6_8h'] <- 1

am_informed['Tin_6_8h', 'Bin_6_8h'] <- 1
am_informed['Bin_6_8h', 'Bin_8_10h'] <- 1
am_informed['Bin_6_8h', 'Bin_10_12h'] <- 1
am_informed['Bin_6_8h', 'Bap_6_8h'] <- 1
am_informed['Bin_8_10h', 'Bin_10_12h'] <- 1
am_informed['Tin_6_8h', 'Bap_6_8h'] <- 1

# Run baycn
# Three independent long runs
baycn_dros <- mhEdge(data = drosophila$continuous,
                     adjMatrix = am_informed,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     nCPh = 0,
                     nGV = 6,
                     pmr = TRUE,
                     burnIn = 0.2,
                     iterations = iterations,
                     thinTo = thinTo,
                     progress = TRUE)

save(baycn_dros, file=paste0("dros_baycn_", iterations, "_", run_number, ".RData"))

cat("Script completed at:", format(Sys.time()), "\n")
