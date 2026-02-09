library (BiDAG)
library (networkBMA)

load('data_geuvadis_pcs.RData')

# adjacency matrix for graphs with nine nodes.
am_9 <- matrix(c(0, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 0, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 0, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 0, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 0, 0, 0, 0, 0,
                 1, 1, 1, 1, 0, 0, 0, 0, 0,
                 1, 1, 1, 1, 0, 0, 0, 0, 0,
                 1, 1, 1, 1, 0, 0, 0, 0, 0,
                 1, 1, 1, 1, 0, 0, 0, 0, 0),
               nrow = 9,
               ncol = 9,
               byrow = TRUE)

blacklist <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = 9,
                    ncol = 9,
                    byrow = TRUE)

# Run order and partition MCMC -------------------------------------------------

# Create the score parameters for order and partition MCMC.
score_q8_pcs <- scoreparameters('bge',
                                  data = data_Q8_pc)

set.seed(1)

# Set MCMC iterations and step sizes
niter.long <- 500000
burn.in <- 0.2
nsamples <- 1000

stepsize.long <- niter.long*(1-burn.in) / nsamples


# Run partition MCMC for the same number of iterations as baycn.
partition_q8_pcs <- partitionMCMC(iterations = niter.long,
                                  startspace = am_9,
                                  scorepar = score_q8_pcs,
                                  blacklist = blacklist,
                                  stepsave = stepsize.long,
                                  verbose = TRUE)

# Calculate the posterior probability adjacency matrix.
partition_posterior <- edgep(partition_q8_pcs,
                             burnin = burn.in)

set.seed(2)

# Run order MCMC for the same number of iterations as baycn.
order_q8_pcs <- orderMCMC(iterations = niter.long,
                          startspace = am_9,
                          scorepar = score_q8_pcs,
                          blacklist = blacklist,
                          stepsave = stepsize.long,
                          chainout = TRUE,
                          MAP = FALSE,
                          verbose = TRUE)

# Calculate the posterior probability adjacency matrix.
order_posterior <- edgep(order_q8_pcs,
                         burnin = burn.in)

# scanBMA helper functions -----------------------------------------------------

# am_bma is a function to create an "adjacency" matrix from the scanBMA output.
am_bma <- function (bma) {
  
  am <- matrix(0,
               nrow = length(bma),
               ncol = length(bma))
  
  # Loop through each node in bma to extract the parent nodes.
  for (e in 1:length(bma)) {
    
    # Create a counter to subset the poterior probability vector from the bma
    # output. This vector is length(bma) - 1 because a node cannot be the parent
    # of itself.
    counter <- 1
    
    # Loop through each potential parent of the current node.
    for (v in 1:length(bma)) {
      
      # Check if the current parent v is the current node e (i.e., the diagonal
      # of the adjacency matrix).
      if (e != v) {
        
        # Add the posterior probability of the current parent to the column of
        # the adjacency matrix corresponding to the current node.
        am[v, e] <- bma[[e]]$probne0[[counter]] / 100
        
        # Increase the counter by one.
        counter <- counter + 1
        
      }
      
    }
    
  }
  
  return (am)
  
}

# Run scanBMA  -------------------------------------------

# The following function will be used on all calls to the scanBMA function
# regardless of which topology is the input.
# Use the default settings (from the networkBMA vignette) for setting the value
# for Occam's window and whether to use Zellner's g prior or BIC.
control <- ScanBMAcontrol(OR = 20,
                          useg = TRUE,
                          gCtrl = gControl(optimize = FALSE,
                                           g0 = 20))

# Create a list that will hold the output from scanBMA for each node.
scan_q8_pcs <- vector(mode = "list",
                      length = dim(am_9)[[1]])

set.seed(3)

# Loop through each variable (scanBMA considers one node/variable at a time).
for (e in 1:dim(am_9)[[1]]) {
  
  scan_q8_pcs[[e]] <- ScanBMA(x = data_Q8_pc[, -e],
                              y = data_Q8_pc[, e],
                              prior.prob = 0.1,
                              control = control)
  
}

scan_posterior <- am_bma(scan_q8_pcs)

# Ignore the posterior probabilities in the first column. We can just zero them
# out because that is the same as if they weren't run. We are doing this because
# we don't want to consider the genetic variant as the child of other nodes.
scan_posterior[2:9, 1] <- 0


# Run BCDAG -----------------------------------------------------
library (BCDAG)
q <- ncol(data_Q8_pc)
n <- nrow(data_Q8_pc)
S_q8 <- 5e5
burn_q8 <- S_q8 * 0.2
thinned_indices<- seq(from=burn_q8+1, to=S_q8, length.out=1000)

bcdag_q8 = learn_DAG(S = S_q8, burn = burn_q8, a = q, U = diag(1,q)/n, data = data_Q8_pc, w = 0.05,
                           fast = FALSE, save.memory = FALSE, collapse = TRUE)
# extract the posterior adjacency matrix
bcdag_thinned <- bcdag_q8$Graphs[,,thinned_indices]
bcdag_q8_adjmtx <- apply(bcdag_thinned, c(1, 2), mean)
colnames(bcdag_q8_adjmtx) <- colnames(data_Q8_pc)
rownames(bcdag_q8_adjmtx) <- colnames(data_Q8_pc)

write.table(bcdag_q8_adjmtx, "bcdag_Q8_adj_mtx.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

# Create edge state matrix -----------------------------------------------------

toPES <- function (pm) {
  
  dims <- dim(pm)
  
  npm <- matrix(as.numeric(pm),
                nrow = dims[1],
                ncol = dims[2])
  
  upper <- t(npm)[lower.tri(npm)]
  lower <- npm[lower.tri(npm)]
  absent <- 1 - (upper + lower)
  
  return (cbind(upper,
                lower,
                absent))
  
}

# The following vectors are the row/column indices of the adjacency matrix that
# correspond to the edges considered by baycn.
the_row <- c(
  1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4
)
the_col <- c(
  2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9
)

order_pes <- matrix(NA, nrow = length(the_row), ncol = 3)
partition_pes <- matrix(NA, nrow = length(the_row), ncol = 3)
scan_pes <- matrix(NA, nrow = length(the_row), ncol = 3)

# Loop through each row/column pair and extract the 0 and 1 edge state values.
for (e in 1:length(the_row)) {
  
  # Extract the value for state 0.
  order_pes[e, 1] <- order_posterior[the_row[[e]], the_col[[e]]]
  partition_pes[e, 1] <- partition_posterior[the_row[[e]], the_col[[e]]]
  scan_pes[e, 1] <- scan_posterior[the_row[[e]], the_col[[e]]]
  
  # Extract the value for state 1.
  order_pes[e, 2] <- order_posterior[the_col[[e]], the_row[[e]]]
  partition_pes[e, 2] <- partition_posterior[the_col[[e]], the_row[[e]]]
  scan_pes[e, 2] <- scan_posterior[the_col[[e]], the_row[[e]]]
  
  # Calculate the value for state 2 (only for order and partition). ScanBMA
  # cannot provide a posterior probability for edge absence.
  order_pes[e, 3] <- 1 - order_pes[e, 1] - order_pes[e, 2]
  partition_pes[e, 3] <- 1 - partition_pes[e, 1] - partition_pes[e, 2]
  
}

# Convert the posterior edge state (PES) matrices to data frames and add the
# node variable.
nodes <- data.frame(
  node = c("rs11305802-TMEM55B", "rs11305802-RP11-203M5.8", "rs11305802-PNP", 
           "rs11305802-PC1", "rs11305802-PC2", "rs11305802-PC6",
           "rs11305802-PC7", "rs11305802-PC9", "TMEM55B-RP11-203M5.8",
           "TMEM55B-PNP", "TMEM55B-PC1", "TMEM55B-PC2", "TMEM55B-PC6",
           "TMEM55B-PC7", "TMEM55B-PC9", "RP11-203M5.8-PNP", "RP11-203M5.8-PC1",
           "RP11-203M5.8-PC2", "RP11-203M5.8-PC6", "RP11-203M5.8-PC7", 
           "RP11-203M5.8-PC9", "PNP-PC1", "PNP-PC2", "PNP-PC6", "PNP-PC7",
           "PNP-PC9")
)

order_pes <- data.frame(
  node = nodes,
  order_pes
)
names(order_pes)[2:4] <- c("zero", "one", "two")

partition_pes <- data.frame(
  node = nodes,
  partition_pes
)
names(partition_pes)[2:4] <- c("zero", "one", "two")

scan_pes <- data.frame(
  node = nodes,
  scan_pes
)
names(scan_pes)[2:4] <- c("zero", "one", "two")

# Save raw output,adjacency and PES matrices -----------------------------------

save(order_q8_pcs,
     partition_q8_pcs,
     scan_q8_pcs,
     
     order_posterior,
     partition_posterior,
     scan_posterior,
     
     order_pes,
     partition_pes,
     scan_pes,
     
     file = "/Users/Evatar/Downloads/ops_geuvadis_q8_pcs.RData")


bcdag_pes <- matrix(NA, nrow = length(the_row), ncol = 3)
# Loop through each row/column pair and extract the 0 and 1 edge state values.
# BCDAG
for (e in 1:length(the_row)) {
  
  # Extract the value for state 0.
  bcdag_pes[e, 1] <- bcdag_q8_adjmtx[the_row[[e]], the_col[[e]]]
  
  # Extract the value for state 1.
  bcdag_pes[e, 2] <- bcdag_q8_adjmtx[the_col[[e]], the_row[[e]]]
  
  # Calculate the value for state 2 (only for order and partition). ScanBMA
  # cannot provide a posterior probability for edge absence.
  bcdag_pes[e, 3] <- 1 - bcdag_pes[e, 1] - bcdag_pes[e, 2]
}

write.table (round (bcdag_pes, digits = 2), "bcdag_q8_pes.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
