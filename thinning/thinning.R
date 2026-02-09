##########################################
# mc3
# no difference with and without thinning
##########################################
# baycn/grid_simulation_analysis
load("mc3_ge_N_600.RData")
tmp <- mc3_g2_600_02_fc[[1]]

tmp <- mc3_g2_600_05_fc[[1]]

# Thin the samples
thinned_indices <- seq(1, 24000, by = 120)
thinned_samples <- tmp$samples[thinned_indices]

# Preserve the class
class(thinned_samples) <- c("mcmcbn", "bn.list", "parental.list")

# Create a modified environment with thinned samples
tmp_thinned <- new.env()
tmp_thinned$samples <- thinned_samples
tmp_thinned$data <- tmp$data
tmp_thinned$logScoreFUN <- tmp$logScoreFUN
tmp_thinned$sampler <- tmp$sampler
tmp_thinned$tabulated <- NULL  # reset this since it's based on original samples
tmp_thinned$type <- tmp$type
class(tmp_thinned) <- "bnpostmcmc"

# Now calculate ep
ep_thinned <- ep(tmp_thinned)
ep_thinned

##########################################
# BCDAG
# Previous runs did not save MCMC samples.
# Rerun with thinning on the grid and save the samples
##########################################
load("bcdag_ge_N_600.RData")

# need to reanalyze the simulated data
setwd("~/Documents/GitHub/baycn_analysis_code/simulation_analysis")
load("data_ge_N_600.RData")

##########################################
# BiDAG (order and partition MCMC)
# Previous orderMCMC runs did not set MAP to FALSE.
# Rerun to get more unbiased samples and allow thinning.
# 
# Previous partitionMCMC runs did not do thinning.
# Rerun to save thinned results.
##########################################
rm(list = ls(pattern = "^mc3_"))
rm(list = ls(pattern = "^out_"))
rm(list = ls(pattern = "^mat_"))

setwd("~/Documents/BreastCancer/baycn/grid_simulation_analysis")

# 100
load("ord_ge_N_100.RData")

# Access the incidence matrices (DAG samples)
samples <- ord_g2_100_02_fc[[1]]$traceadd$incidence
length(samples)  # 1001
# This is because the default for the stepsave argument 
# in orderMCMC() is 1000

# Look at one sample
samples[[1]]

# Thin to ~200 samples (every 5th)
thinned_indices <- seq(1, length(samples), by = 5)
thinned_samples <- samples[thinned_indices]
length(thinned_samples)  # ~200

# Calculate posterior edge probabilities manually
# Convert list of matrices to 3D array and take mean
# Convert all to regular matrices
thinned_matrices <- lapply(thinned_samples, as.matrix)

# Now create the 3D array and calculate mean
sample_array <- array(unlist(thinned_matrices), 
                      dim = c(4, 4, length(thinned_matrices)))
ep_thinned <- apply(sample_array, c(1, 2), mean)
ep_thinned


calc_ep_bidag <- function(mcmc_list, burnin = 0, thin = 5) {
  # Calculate posterior edge probabilities from BiDAG orderMCMC output
  #
  # Args:
  #   mcmc_list: list of orderMCMC output objects
  #   burnin: proportion of samples to discard as burn-in (default 0.2)
  #   thin: thinning interval (default 5, take every 5th sample)
  #
  # Returns:
  #   list of posterior adjacency matrices (same length as input)
  
  n_chains <- length(mcmc_list)
  ep_list <- vector("list", n_chains)
  
  for (i in 1:n_chains) {
    # Extract incidence matrices
    samples <- mcmc_list[[i]]$traceadd$incidence
    n_samples <- length(samples)
    
    # Get matrix dimension from first sample
    p <- nrow(as.matrix(samples[[1]]))
    
    # Apply burn-in
    burnin_n <- floor(burnin * n_samples)
    post_burnin <- samples[(burnin_n + 1):n_samples]
    
    # Thin the samples
    thinned_indices <- seq(1, length(post_burnin), by = thin)
    thinned_samples <- post_burnin[thinned_indices]
    
    # Convert all to regular matrices (handles mix of sparse and dense)
    thinned_matrices <- lapply(thinned_samples, as.matrix)
    
    # Create 3D array and calculate mean
    sample_array <- array(unlist(thinned_matrices), 
                          dim = c(p, p, length(thinned_matrices)))
    ep_list[[i]] <- apply(sample_array, c(1, 2), mean)
  }
  
  return(ep_list)
}

# Process bcdag output
files <- list.files("./bcdag_output", pattern = ".*g2_600_02.*\\.rds$", full.names = TRUE)
bcdag_g2_600_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*g2_600_05.*\\.rds$", full.names = TRUE)
bcdag_g2_600_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*g2_600_1.*\\.rds$", full.names = TRUE)
bcdag_g2_600_1_pm <- lapply(files, readRDS)

files <- list.files("./bcdag_output", pattern = ".*nc11_600_02.*\\.rds$", full.names = TRUE)
bcdag_nc11_600_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*nc11_600_05.*\\.rds$", full.names = TRUE)
bcdag_nc11_600_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*nc11_600_1.*\\.rds$", full.names = TRUE)
bcdag_nc11_600_1_pm <- lapply(files, readRDS)

files <- list.files("./bcdag_output", pattern = ".*pc_600_02.*\\.rds$", full.names = TRUE)
bcdag_pc_600_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*pc_600_05.*\\.rds$", full.names = TRUE)
bcdag_pc_600_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*pc_600_1.*\\.rds$", full.names = TRUE)
bcdag_pc_600_1_pm <- lapply(files, readRDS)

save(bcdag_g2_600_02_pm,
     bcdag_g2_600_05_pm,
     bcdag_g2_600_1_pm,
     bcdag_nc11_600_02_pm,
     bcdag_nc11_600_05_pm,
     bcdag_nc11_600_1_pm,
     bcdag_pc_600_02_pm,
     bcdag_pc_600_05_pm,
     bcdag_pc_600_1_pm,
     file = "./bcdag_ge_N_600_pm.RData"
)

files <- list.files("./bcdag_output", pattern = ".*g2_200_02.*\\.rds$", full.names = TRUE)
bcdag_g2_200_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*g2_200_05.*\\.rds$", full.names = TRUE)
bcdag_g2_200_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*g2_200_1.*\\.rds$", full.names = TRUE)
bcdag_g2_200_1_pm <- lapply(files, readRDS)

files <- list.files("./bcdag_output", pattern = ".*nc11_200_02.*\\.rds$", full.names = TRUE)
bcdag_nc11_200_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*nc11_200_05.*\\.rds$", full.names = TRUE)
bcdag_nc11_200_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*nc11_200_1.*\\.rds$", full.names = TRUE)
bcdag_nc11_200_1_pm <- lapply(files, readRDS)

files <- list.files("./bcdag_output", pattern = ".*pc_200_02.*\\.rds$", full.names = TRUE)
bcdag_pc_200_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*pc_200_05.*\\.rds$", full.names = TRUE)
bcdag_pc_200_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*pc_200_1.*\\.rds$", full.names = TRUE)
bcdag_pc_200_1_pm <- lapply(files, readRDS)

save(bcdag_g2_200_02_pm,
     bcdag_g2_200_05_pm,
     bcdag_g2_200_1_pm,
     bcdag_nc11_200_02_pm,
     bcdag_nc11_200_05_pm,
     bcdag_nc11_200_1_pm,
     bcdag_pc_200_02_pm,
     bcdag_pc_200_05_pm,
     bcdag_pc_200_1_pm,
     file = "./bcdag_ge_N_200_pm.RData"
)

# 100
files <- list.files("./bcdag_output", pattern = ".*g2_100_02.*\\.rds$", full.names = TRUE)
bcdag_g2_100_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*g2_100_05.*\\.rds$", full.names = TRUE)
bcdag_g2_100_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*g2_100_1.*\\.rds$", full.names = TRUE)
bcdag_g2_100_1_pm <- lapply(files, readRDS)

files <- list.files("./bcdag_output", pattern = ".*nc11_100_02.*\\.rds$", full.names = TRUE)
bcdag_nc11_100_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*nc11_100_05.*\\.rds$", full.names = TRUE)
bcdag_nc11_100_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*nc11_100_1.*\\.rds$", full.names = TRUE)
bcdag_nc11_100_1_pm <- lapply(files, readRDS)

files <- list.files("./bcdag_output", pattern = ".*pc_100_02.*\\.rds$", full.names = TRUE)
bcdag_pc_100_02_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*pc_100_05.*\\.rds$", full.names = TRUE)
bcdag_pc_100_05_pm <- lapply(files, readRDS)
files <- list.files("./bcdag_output", pattern = ".*pc_100_1.*\\.rds$", full.names = TRUE)
bcdag_pc_100_1_pm <- lapply(files, readRDS)

save(bcdag_g2_100_02_pm,
     bcdag_g2_100_05_pm,
     bcdag_g2_100_1_pm,
     bcdag_nc11_100_02_pm,
     bcdag_nc11_100_05_pm,
     bcdag_nc11_100_1_pm,
     bcdag_pc_100_02_pm,
     bcdag_pc_100_05_pm,
     bcdag_pc_100_1_pm,
     file = "./bcdag_ge_N_100_pm.RData"
)

# Calculate MSE2, precision and recall
setwd("~/Documents/BreastCancer/baycn/thinning")
load("./ord_ge_N_600_pm.RData")
load("./part_ge_N_600_pm.RData")
load("./bcdag_ge_N_600_pm.RData")

# expected probability matrices ------------------------------------------------

# The expected probabilities for topology G2.
ep_g2 <- matrix(c(0, 1/3, 1, 0,
                  2/3, 0, 0, 2/3,
                  0, 0, 0, 0,
                  0, 1/3, 1, 0),
                byrow = TRUE,
                nrow = 4)

# The expected probabilities for topology NC11.
ep_nc11 <- matrix(c(0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0.8, 0, 0.4, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0.6, 0, 0.6, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0.4, 0, 0.8, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0.2, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 1, 0, 0.2, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0.8, 0, 0.4, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0.6, 0, 0.6, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0.4, 0, 0.8,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2, 0),
                  byrow = TRUE,
                  nrow = 11)

# The expected probabilities for topology PC.
ep_pc <- matrix(c(0, 0.25, 0, 0, 0, 1, 0, 1,
                  0.75, 0, 0.75, 0, 0.75, 0, 0, 0,
                  0, 0.25, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0.25, 0, 0, 0, 1, 0, 1,
                  0, 0, 0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0),
                byrow = TRUE,
                nrow = 8)

# The MSE is calculated by subtracting the estimated edge probabilities (not
# including the probability of the edge being absent) from the expected edge
# probabilities.

calc_mse <- function(pm_list, true_mat) {
  n_nodes <- nrow (true_mat)
  n_edges <- n_nodes * (n_nodes - 1)  # off-diagonal elements
  sapply(pm_list, function(x) sum((true_mat - x)^2) / n_edges)
}

ord_mse_g2_600_02 <- calc_mse(ord_g2_600_02_pm, ep_g2)
ord_mse_g2_600_05 <- calc_mse(ord_g2_600_05_pm, ep_g2)
ord_mse_g2_600_1 <- calc_mse(ord_g2_600_1_pm, ep_g2)

ord_mse_nc11_600_02 <- calc_mse(ord_nc11_600_02_pm, ep_nc11)
ord_mse_nc11_600_05 <- calc_mse(ord_nc11_600_05_pm, ep_nc11)
ord_mse_nc11_600_1 <- calc_mse(ord_nc11_600_1_pm, ep_nc11)

ord_mse_pc_600_02 <- calc_mse(ord_pc_600_02_pm, ep_pc)
ord_mse_pc_600_05 <- calc_mse(ord_pc_600_05_pm, ep_pc)
ord_mse_pc_600_1 <- calc_mse(ord_pc_600_1_pm, ep_pc)

part_mse_g2_600_02 <- calc_mse(part_g2_600_02_pm, ep_g2)
part_mse_g2_600_05 <- calc_mse(part_g2_600_05_pm, ep_g2)
part_mse_g2_600_1 <- calc_mse(part_g2_600_1_pm, ep_g2)

part_mse_nc11_600_02 <- calc_mse(part_nc11_600_02_pm, ep_nc11)
part_mse_nc11_600_05 <- calc_mse(part_nc11_600_05_pm, ep_nc11)
part_mse_nc11_600_1 <- calc_mse(part_nc11_600_1_pm, ep_nc11)

part_mse_pc_600_02 <- calc_mse(part_pc_600_02_pm, ep_pc)
part_mse_pc_600_05 <- calc_mse(part_pc_600_05_pm, ep_pc)
part_mse_pc_600_1 <- calc_mse(part_pc_600_1_pm, ep_pc)


mean (unlist (part_mse_g2_600_02))
sd (unlist (part_mse_g2_600_02))
mean (unlist (part_mse_g2_600_05))
sd (unlist (part_mse_g2_600_05))
mean (unlist (part_mse_g2_600_1))
sd (unlist (part_mse_g2_600_1))

mean (unlist (part_mse_nc11_600_02))
sd (unlist (part_mse_nc11_600_02))
mean (unlist (part_mse_nc11_600_05))
sd (unlist (part_mse_nc11_600_05))
mean (unlist (part_mse_nc11_600_1))
sd (unlist (part_mse_nc11_600_1))

mean (unlist (part_mse_pc_600_02))
sd (unlist (part_mse_pc_600_02))
mean (unlist (part_mse_pc_600_05))
sd (unlist (part_mse_pc_600_05))
mean (unlist (part_mse_pc_600_1))
sd (unlist (part_mse_pc_600_1))


bcdag_mse_g2_600_02 <- calc_mse(bcdag_g2_600_02_pm, ep_g2)
bcdag_mse_g2_600_05 <- calc_mse(bcdag_g2_600_05_pm, ep_g2)
bcdag_mse_g2_600_1 <- calc_mse(bcdag_g2_600_1_pm, ep_g2)

bcdag_mse_nc11_600_02 <- calc_mse(bcdag_nc11_600_02_pm, ep_nc11)
bcdag_mse_nc11_600_05 <- calc_mse(bcdag_nc11_600_05_pm, ep_nc11)
bcdag_mse_nc11_600_1 <- calc_mse(bcdag_nc11_600_1_pm, ep_nc11)

bcdag_mse_pc_600_02 <- calc_mse(bcdag_pc_600_02_pm, ep_pc)
bcdag_mse_pc_600_05 <- calc_mse(bcdag_pc_600_05_pm, ep_pc)
bcdag_mse_pc_600_1 <- calc_mse(bcdag_pc_600_1_pm, ep_pc)

mean (unlist (bcdag_mse_g2_600_02))
sd (unlist (bcdag_mse_g2_600_02))
mean (unlist (bcdag_mse_g2_600_05))
sd (unlist (bcdag_mse_g2_600_05))
mean (unlist (bcdag_mse_g2_600_1))
sd (unlist (bcdag_mse_g2_600_1))

mean (unlist (bcdag_mse_nc11_600_02))
sd (unlist (bcdag_mse_nc11_600_02))
mean (unlist (bcdag_mse_nc11_600_05))
sd (unlist (bcdag_mse_nc11_600_05))
mean (unlist (bcdag_mse_nc11_600_1))
sd (unlist (bcdag_mse_nc11_600_1))

mean (unlist (bcdag_mse_pc_600_02))
sd (unlist (bcdag_mse_pc_600_02))
mean (unlist (bcdag_mse_pc_600_05))
sd (unlist (bcdag_mse_pc_600_05))
mean (unlist (bcdag_mse_pc_600_1))
sd (unlist (bcdag_mse_pc_600_1))
