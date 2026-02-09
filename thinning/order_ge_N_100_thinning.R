cat("Script started at:", format(Sys.time()), "\n")

# Define home directory
home_dir <- "/wsu/home/ht/ht26/ht2699"

# Load libraries
library ("BiDAG", lib=file.path(home_dir, "Rpackages"))
#library ("gRbase", lib=file.path(home_dir, "Rpackages"))
#library ("mvtnorm", lib=file.path(home_dir, "Rpackages"))

# Load the data
load(file.path(home_dir, "baycn_data_analysis/simulation/data/data_ge_N_100.RData"))

# Set M value
M <- 25

# Set up the starting adjacency matrices
# Undirected adjacency matrix with the all edges for topology G2.
am_g2 <- matrix(c(0, 1, 1, 1,
                  1, 0, 1, 1,
                  1, 1, 0, 1,
                  1, 1, 1, 0),
                byrow = TRUE,
                nrow = 4)

# Undirected adjacency matrix with the all edges for toplogy NC11.
am_nc11 <- matrix(c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0),
                  byrow = TRUE,
                  nrow = 11)

# Undirected adjacency matrix with the all edges for the PC topology.
am_pc <- matrix(c(0, 1, 1, 1, 1, 1, 1, 1,
                  1, 0, 1, 1, 1, 1, 1, 1,
                  1, 1, 0, 1, 1, 1, 1, 1,
                  1, 1, 1, 0, 1, 1, 1, 1,
                  1, 1, 1, 1, 0, 1, 1, 1,
                  1, 1, 1, 1, 1, 0, 1, 1,
                  1, 1, 1, 1, 1, 1, 0, 1,
                  1, 1, 1, 1, 1, 1, 1, 0),
                byrow = TRUE,
                nrow = 8)


# Initialize the lists for the probability matrix to full length.
ord_g2_100_02_pm <- vector(mode = 'list',
                           length = M)

ord_g2_100_05_pm <- vector(mode = 'list',
                           length = M)

ord_g2_100_1_pm <- vector(mode = 'list',
                          length = M)

ord_nc11_100_02_pm <- vector(mode = 'list',
                             length = M)

ord_nc11_100_05_pm <- vector(mode = 'list',
                             length = M)

ord_nc11_100_1_pm <- vector(mode = 'list',
                            length = M)

ord_pc_100_02_pm <- vector(mode = 'list',
                           length = M)

ord_pc_100_05_pm <- vector(mode = 'list',
                           length = M)

ord_pc_100_1_pm <- vector(mode = 'list',
                          length = M)

# Initialize the time lists to full length.
ord_g2_100_02_time <- vector(mode = 'list',
                             length = M)

ord_g2_100_05_time <- vector(mode = 'list',
                             length = M)

ord_g2_100_1_time <- vector(mode = 'list',
                            length = M)

ord_nc11_100_02_time <- vector(mode = 'list',
                               length = M)

ord_nc11_100_05_time <- vector(mode = 'list',
                               length = M)

ord_nc11_100_1_time <- vector(mode = 'list',
                              length = M)

ord_pc_100_02_time <- vector(mode = 'list',
                             length = M)

ord_pc_100_05_time <- vector(mode = 'list',
                             length = M)

ord_pc_100_1_time <- vector(mode = 'list',
                            length = M)

#set.seed(200)

# Set MCMC iterations and step sizes
niter.short <- 30000
niter.long <- 50000

burn.in <- 0.2
nsamples <- 200

stepsize.short <- niter.short*(1-burn.in) / nsamples
stepsize.long <- niter.long*(1-burn.in) / nsamples

# Loop through each combination of signal strength and sample size for all
# topologies.
for (e in 1:M) {

  #################################################
  # G2
  # 100
  # 0.2
  #################################################
  cat("Processing G2, replicate", e, "of", M, "\n")
  
  score_g2_100_02 <- scoreparameters("bge",
                                     data_g2_100_02[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.short,
                                  scorepar = score_g2_100_02,
                                  startspace = am_g2,
                                  stepsave = stepsize.short, 
                                  chainout = TRUE,
                                  MAP = FALSE)
  
  endtime <- Sys.time()
  ord_g2_100_02_time[[e]] <- endtime - starttime
  ord_g2_100_02_pm[[e]] <- edgep(temp_out, burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # G2 
  # 100
  # 0.5
  #################################################

  score_g2_100_05 <- scoreparameters(
                                     "bge",
                                     data_g2_100_05[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.short,
                                  scorepar = score_g2_100_05,
                                  startspace = am_g2,
                                  stepsave = stepsize.short, 
                                  chainout = TRUE,
                                  MAP = FALSE)
  endtime <- Sys.time()
  ord_g2_100_05_time[[e]] <- endtime - starttime
  ord_g2_100_05_pm[[e]] <- edgep(temp_out,
                                           burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # G2 
  # 100
  # 1
  #################################################

  score_g2_100_1 <- scoreparameters(
                                    "bge",
                                    data_g2_100_1[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.short,
                                  scorepar = score_g2_100_1,
                                  startspace = am_g2,
                                  stepsave = stepsize.short, 
                                  chainout = TRUE,
                                  MAP = FALSE)
  endtime <- Sys.time()
  ord_g2_100_1_time[[e]] <- endtime - starttime
  ord_g2_100_1_pm[[e]] <- edgep(temp_out,
                                          burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # NC11
  # 100
  # 0.2
  #################################################

  cat("Processing NC11, replicate", e, "of", M, "\n")

  score_nc11_100_02 <- scoreparameters(
                                       "bge",
                                       data_nc11_100_02[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.long,
                                  scorepar = score_nc11_100_02,
                                  startspace = am_nc11,
                                  stepsave = stepsize.long, 
                                  chainout = TRUE,
                                  MAP = FALSE)
  endtime <- Sys.time()
  ord_nc11_100_02_time[[e]] <- endtime - starttime
  ord_nc11_100_02_pm[[e]] <- edgep(temp_out,
                                             burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # NC11
  # 100
  # 0.5
  #################################################

  score_nc11_100_05 <- scoreparameters(
                                       "bge",
                                       data_nc11_100_05[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.long,
                                    scorepar = score_nc11_100_05,
                                    startspace = am_nc11,
                                    stepsave = stepsize.long, 
                                    chainout = TRUE,
                                    MAP = FALSE)
  endtime <- Sys.time()
  ord_nc11_100_05_time[[e]] <- endtime - starttime
  ord_nc11_100_05_pm[[e]] <- edgep(temp_out,
                                             burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # NC11
  # 100
  # 1
  #################################################

  score_nc11_100_1 <- scoreparameters(
                                      "bge",
                                      data_nc11_100_1[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.long,
                                    scorepar = score_nc11_100_1,
                                    startspace = am_nc11,
                                    stepsave = stepsize.long, 
                                    chainout = TRUE,
                                    MAP = FALSE)
  endtime <- Sys.time()
  ord_nc11_100_1_time[[e]] <- endtime - starttime
  ord_nc11_100_1_pm[[e]] <- edgep(temp_out,
                                             burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # PC
  # 100
  # 0.2
  #################################################
  cat("Processing PC, replicate", e, "of", M, "\n")
  
  score_pc_100_02 <- scoreparameters(
                                     "bge",
                                     data_pc_100_02[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.long,
                                    scorepar = score_pc_100_02,
                                    startspace = am_pc,
                                    stepsave = stepsize.long, 
                                    chainout = TRUE,
                                    MAP = FALSE)
  endtime <- Sys.time()
  ord_pc_100_02_time[[e]] <- endtime - starttime
  ord_pc_100_02_pm[[e]] <- edgep(temp_out,
                                           burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # PC
  # 100
  # 0.5
  #################################################

  score_pc_100_05 <- scoreparameters(
                                     "bge",
                                     data_pc_100_05[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.long,
                                  scorepar = score_pc_100_05,
                                  startspace = am_pc,
                                  stepsave = stepsize.long, 
                                  chainout = TRUE,
                                  MAP = FALSE)
  endtime <- Sys.time()
  ord_pc_100_05_time[[e]] <- endtime - starttime
  ord_pc_100_05_pm[[e]] <- edgep(temp_out,
                                           burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)
  
  #################################################
  # PC
  # 100
  # 1
  #################################################

  score_pc_100_1 <- scoreparameters(
                                    "bge",
                                    data_pc_100_1[[e]])

  starttime <- Sys.time()
  temp_out <- orderMCMC(iterations = niter.long,
                                  scorepar = score_pc_100_1,
                                  startspace = am_pc,
                                  stepsave = stepsize.long, 
                                  chainout = TRUE,
                                  MAP = FALSE)
  endtime <- Sys.time()
  ord_pc_100_1_time[[e]] <- endtime - starttime
  ord_pc_100_1_pm[[e]] <- edgep(temp_out,
                                           burnin = burn.in)
  rm(temp_out); gc(verbose=FALSE)

  cat("Completed replicate", e, "at:", format(Sys.time()), "\n")

  # Save after each replicate for crash recovery
  save(M,
       ord_g2_100_02_pm,
       ord_g2_100_05_pm,
       ord_g2_100_1_pm,
       ord_nc11_100_02_pm,
       ord_nc11_100_05_pm,
       ord_nc11_100_1_pm,
       ord_pc_100_02_pm,
       ord_pc_100_05_pm,
       ord_pc_100_1_pm,
       file = './ord_ge_N_100_pm.RData')

  save(M,
       ord_g2_100_02_time,
       ord_g2_100_05_time,
       ord_g2_100_1_time,
       ord_nc11_100_02_time,
       ord_nc11_100_05_time,
       ord_nc11_100_1_time,
       ord_pc_100_02_time,
       ord_pc_100_05_time,
       ord_pc_100_1_time,
       file = './ord_ge_N_100_time.RData')

}

cat("Script completed at:", format(Sys.time()), "\n")
