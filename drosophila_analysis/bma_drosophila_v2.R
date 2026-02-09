library("RcppArmadillo", lib="~/bin/r.packages/3.2.3")
library ("networkBMA", lib="~/bin/r.packages/3.2.3")

# Load the drosphila data from baycn ----------------------------------
# data("drosophila")
# load the saved drosophila from .csv
drosophila.continuous = read.csv("./drosophila.continuous.csv")

# Start time 
start_time <- Sys.time()

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
scan_continuous <- vector(mode = "list",
                        length = 21) # changed from dim(am_symmetric)[[1]] to 21

set.seed(19)

# Loop through each variable (scanBMA considers one node/variable at a time).
for (e in 1:21) {
  
  scan_continuous[[e]] <- ScanBMA(x = drosophila.continuous[, -e],
                                y = drosophila.continuous[, e],
                                prior.prob = 0.1,
                                control = control)
  
}

scan_posterior <- am_bma(scan_continuous)

end_time = Sys.time()
total_time = end_time - start_time

print(total_time)

save(scan_posterior, file="drosophila.continuous.scan_posterior.RData")
save(total_time, file="drosophila.continuous.scan_posterior.runtime.RData")
print("Complete ScanBMA run")
