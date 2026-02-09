# Load necessary packages ------------------------------------------------------

library (BiDAG)
library (networkBMA)

# Load the drosphila data from baycn ----------------------------------
data("drosophila")
drosophila_discrete <- drosophila$discrete
drosophila_continuous <- drosophila$continuous

# Load adjacency matrix with additional edges ----------------------------------

# Load matrix from the same directory as this script
load('am_discrete.RData')

# Run order and partition MCMC on the drosophila data --------------------------

# Create the search space for order and partition MCMC. This has to be a
# symmetric matrix in order for both directions of each edge to be considered.
am_symmetric <- ifelse((am_discrete + t(am_discrete)) > 0, 1, 0)

# Create the score parameters for order and partition MCMC.
score_discrete <- scoreparameters('bde',
                                  data = drosophila_discrete)

set.seed(21)

# Run partition MCMC for the same number of iterations as baycn.
partition_discrete <- partitionMCMC(iterations = 500000,
                                    startspace = am_symmetric,
                                    scorepar = score_discrete,
                                    verbose = TRUE)

# Calculate the posterior probability adjacency matrix.
partition_posterior <- edgep(partition_discrete,
                             burnin = 0.2)

# Convert the posterior probability adjacency matrix into an edge state matrix.
# partition_pes <- toPES(partition_posterior, am_partition)

# Save the partition MCMC output, posterior probability adjacency matrix, and
# posterior edge state matrix.
# save(partition_discrete,
#      partition_posterior,
#      partition_pes,
#      file = '/Users/Evatar/Desktop/par_drosophila.RData')

set.seed(20)

# Run order MCMC for the same number of iterations as baycn.
order_discrete <- orderMCMC(iterations = 500000,
                            startspace = am_symmetric,
                            scorepar = score_discrete,
                            chainout = TRUE,
                            verbose = TRUE)

# Calculate the posterior probability adjacency matrix.
order_posterior <- edgep(order_discrete,
                         burnin = 0.2)

# Convert the posterior probability adjacency matrix into an edge state matrix.
# order_pes <- toPES(order_posterior, am_symmetric)

# Save the order MCMC output, posterior probability adjacency matrix, and
# posterior edge state matrix.
# save(order_discrete,
#      order_posterior,
#      order_pes,
#      file = '/Users/Evatar/Documents/ord_drosophila.RData')


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

# Run scanBMA on the drosophila data -------------------------------------------

# The following function will be used on all calls to the scanBMA function
# regardless of which topology is the input.
# Use the default settings (from the networkBMA vignette) for setting the value
# for Occam's window and whether to use Zellner's g prior or BIC.
control <- ScanBMAcontrol(OR = 20,
                          useg = TRUE,
                          gCtrl = gControl(optimize = FALSE,
                                           g0 = 20))

# Create a list that will hold the output from scanBMA for each node.
scan_discrete <- vector(mode = "list",
                        length = dim(am_symmetric)[[1]])

set.seed(19)

# Loop through each variable (scanBMA considers one node/variable at a time).
for (e in 1:dim(am_discrete)[[1]]) {
  
  scan_discrete[[e]] <- ScanBMA(x = drosophila$discrete[, -e],
                                y = drosophila$discrete[, e],
                                prior.prob = 0.1,
                                control = control)
  
}

scan_posterior <- am_bma(scan_discrete)


# BCDAG
# BCDAG requires continuous data
# and is a stretch for the mixed drosophila data
library (BCDAG)
q <- ncol(drosophila$continuous)
n <- nrow(drosophila$continuous)
bcdag_discrete = learn_DAG(S = 50000, burn = 10000, a = q, U = diag(1,q)/n, data = drosophila$continuous, w = 0.05,
                           fast = FALSE, save.memory = FALSE, collapse = TRUE)
# extract the posterior adjacency matrix
bcdag_discrete_adjmtx <- get_edgeprobs(bcdag_discrete)
colnames(bcdag_discrete_adjmtx) <- colnames(drosophila$continuous)
rownames(bcdag_discrete_adjmtx) <- colnames(drosophila$continuous)

write.table(bcdag_discrete_adjmtx, "bcdag_drosophila_adj_mtx.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

# Save raw output and adjacency matrices ---------------------------------------

save(order_discrete,
     partition_discrete,
     scan_discrete,
     
     order_posterior,
     partition_posterior,
     scan_posterior,
     
     file = "/Users/Evatar/Documents/ops_drosophila.RData")

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
  1, 2, 2, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8,
  8, 9, 9, 10, 10, 10, 11, 13, 14, 14, 15, 17, 17, 18, 19, 20
)
the_col <- c(
  13, 10, 11, 12, 21, 12, 7, 8, 9, 17, 18, 19, 20, 21, 12, 13, 21, 8, 14,
  14, 15, 17, 18, 10, 16, 11, 12, 16, 12, 14, 15, 18, 19, 19, 21, 19, 20, 21
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

bcdag_pes <- matrix(NA, nrow = length(the_row), ncol = 3)
# Loop through each row/column pair and extract the 0 and 1 edge state values.
# BCDAG
for (e in 1:length(the_row)) {
  
  # Extract the value for state 0.
  bcdag_pes[e, 1] <- bcdag_discrete_adjmtx[the_row[[e]], the_col[[e]]]

  # Extract the value for state 1.
  bcdag_pes[e, 2] <- bcdag_discrete_adjmtx[the_col[[e]], the_row[[e]]]

  # Calculate the value for state 2 (only for order and partition). ScanBMA
  # cannot provide a posterior probability for edge absence.
  bcdag_pes[e, 3] <- 1 - bcdag_pes[e, 1] - bcdag_pes[e, 2]
}
