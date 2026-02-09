cat("Script started at:", format(Sys.time()), "\n")

# Define home directory
home_dir <- "/wsu/home/ht/ht26/ht2699"

# Load libraries
library ("BCDAG", lib=file.path(home_dir, "Rpackages"))
library ("gRbase", lib=file.path(home_dir, "Rpackages"))
library ("mvtnorm", lib=file.path(home_dir, "Rpackages"))

# Load the data
load(file.path(home_dir, "baycn_data_analysis/simulation/data/data_ge_N_100.RData"))

# Set M value
M <- 25

# Initialize model parameters
# a is the shape parameter
w = 0.05
n = 100
a_g2 = 4
a_nc11 = 11 
a_pc = 8
U_g2 = diag(1,a_g2)/n
U_nc11 = diag(1,a_nc11)/n
U_pc = diag(1,a_pc)/n

# Initialize S and burn for the three graphs: g2, nc11, pc
S_g2 <- 30000
burn_g2 <- S_g2*0.2

S_nc11 <- 50000
burn_nc11 <- S_nc11*0.2

S_pc <- 50000
burn_pc <- S_pc*0.2

thinned_indices_short <- seq(from=burn_g2+1, to=S_g2, length.out=200)
thinned_indices_long <- seq(from=burn_nc11+1, to=S_nc11, length.out=200)

# Initiate lists ONLY for matrices (not full outputs) to save memory
# mat_ will store the edge probabilities
mat_g2_100_02 <- vector(mode="list", length=M)
mat_g2_100_05 <- vector(mode="list", length=M)
mat_g2_100_1 <- vector(mode="list", length=M)

mat_nc11_100_02 <- vector(mode="list", length=M)
mat_nc11_100_05 <- vector(mode="list", length=M)
mat_nc11_100_1 <- vector(mode="list", length=M)

mat_pc_100_02 <- vector(mode="list", length=M)
mat_pc_100_05 <- vector(mode="list", length=M)
mat_pc_100_1 <- vector(mode="list", length=M)

# Time tracking
out_g2_100_02_time <- vector(mode="list", length=M)
out_g2_100_05_time <- vector(mode="list", length=M)
out_g2_100_1_time <- vector(mode="list", length=M)

out_nc11_100_02_time <- vector(mode="list", length=M)
out_nc11_100_05_time <- vector(mode="list", length=M)
out_nc11_100_1_time <- vector(mode="list", length=M)

out_pc_100_02_time <- vector(mode="list", length=M)
out_pc_100_05_time <- vector(mode="list", length=M)
out_pc_100_1_time <- vector(mode="list", length=M)


# Loop through each combination of signal strength and sample size for all topologies
for (e in 1:M) {

  ####################################
  # G2
  ####################################

  cat("Processing G2, replicate", e, "of", M, "\n")

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_g2, burn=burn_g2, data=data_g2_100_02[[e]], a=a_g2, U=U_g2, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_short]
  mat_g2_100_02[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_g2_100_02_time[[e]] = endtime - starttime
  saveRDS(mat_g2_100_02[[e]], file=sprintf("./mat_g2_100_02_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_g2, burn=burn_g2, data=data_g2_100_05[[e]], a=a_g2, U=U_g2, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_short]
  mat_g2_100_05[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_g2_100_05_time[[e]] = endtime - starttime
  saveRDS(mat_g2_100_05[[e]], file=sprintf("./mat_g2_100_05_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_g2, burn=burn_g2, data=data_g2_100_1[[e]], a=a_g2, U=U_g2, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_short]
  mat_g2_100_1[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_g2_100_1_time[[e]] = endtime - starttime
  saveRDS(mat_g2_100_1[[e]], file=sprintf("./mat_g2_100_1_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  ####################################
  # NC11
  ####################################

  cat("Processing NC11, replicate", e, "of", M, "\n")

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_nc11, burn=burn_nc11, data=data_nc11_100_02[[e]], a=a_nc11, U=U_nc11, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_long]
  mat_nc11_100_02[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_nc11_100_02_time[[e]] = endtime - starttime
  saveRDS(mat_nc11_100_02[[e]], file=sprintf("./mat_nc11_100_02_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_nc11, burn=burn_nc11, data=data_nc11_100_05[[e]], a=a_nc11, U=U_nc11, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_long]
  mat_nc11_100_05[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_nc11_100_05_time[[e]] = endtime - starttime
  saveRDS(mat_nc11_100_05[[e]], file=sprintf("./mat_nc11_100_05_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_nc11, burn=burn_nc11, data=data_nc11_100_1[[e]], a=a_nc11, U=U_nc11, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_long]
  mat_nc11_100_1[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_nc11_100_1_time[[e]] = endtime - starttime
  saveRDS(mat_nc11_100_1[[e]], file=sprintf("./mat_nc11_100_1_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  ####################################
  # PC
  ####################################

  cat("Processing PC, replicate", e, "of", M, "\n")

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_pc, burn=burn_pc, data=data_pc_100_02[[e]], a=a_pc, U=U_pc, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_long]
  mat_pc_100_02[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_pc_100_02_time[[e]] = endtime - starttime
  saveRDS(mat_pc_100_02[[e]], file=sprintf("./mat_pc_100_02_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_pc, burn=burn_pc, data=data_pc_100_05[[e]], a=a_pc, U=U_pc, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_long]
  mat_pc_100_05[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_pc_100_05_time[[e]] = endtime - starttime
  saveRDS(mat_pc_100_05[[e]], file=sprintf("./mat_pc_100_05_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

  starttime <- Sys.time()
  temp_out = learn_DAG(S=S_pc, burn=burn_pc, data=data_pc_100_1[[e]], a=a_pc, U=U_pc, w=w, verbose=FALSE)
  endtime <- Sys.time()
  bcdag_thinned <- temp_out$Graphs[,,thinned_indices_long]
  mat_pc_100_1[[e]] <- apply(bcdag_thinned, c(1, 2), mean)
  out_pc_100_1_time[[e]] = endtime - starttime
  saveRDS(mat_pc_100_1[[e]], file=sprintf("./mat_pc_100_1_rep%02d.rds", e))
  rm(temp_out, bcdag_thinned); gc(verbose=FALSE)

}

# Final save - combine all results into one RData file
cat("All replicates completed. Final save...\n")
save(M,
     mat_g2_100_02, mat_g2_100_05, mat_g2_100_1,
     mat_nc11_100_02, mat_nc11_100_05, mat_nc11_100_1,
     mat_pc_100_02, mat_pc_100_05, mat_pc_100_1,
     out_g2_100_02_time, out_g2_100_05_time, out_g2_100_1_time,
     out_nc11_100_02_time, out_nc11_100_05_time, out_nc11_100_1_time,
     out_pc_100_02_time, out_pc_100_05_time, out_pc_100_1_time,
     file = "./bcdag_ge_N_100_thin.RData"
)

cat("Script completed at:", format(Sys.time()), "\n")

