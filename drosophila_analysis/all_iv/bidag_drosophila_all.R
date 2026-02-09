# Load necessary packages ------------------------------------------------------

library (BiDAG)

# Adjacency matrices 
am_all <- matrix(1, nrow = 21, ncol = 21)
diag(am_all) <- 0
am_all[,1:6] <- 0
rownames(am_all) <- colnames(drosophila$continuous)
colnames(am_all) <- rownames(am_all)

# Generate blacklist to prevent edges from TF to tissue type
blacklist <- matrix(0, nrow = 21, ncol = 21)
rownames(blacklist) <- colnames(drosophila$continuous)
colnames(blacklist) <- rownames(blacklist)
blacklist[,1:6] <- 1

# Set MCMC iterations and step sizes
niter.long <- 500000

burn.in <- 0.2
nsamples <- 1000

stepsize.long <- niter.long*(1-burn.in) / nsamples

# Partition MCMC
score_dros <- scoreparameters("bge", drosophila$continuous)

temp_out <- partitionMCMC(iterations = niter.long,
                          scorepar = score_dros,
                          startspace = am_all,
                          blacklist = blacklist,
                          stepsave = stepsize.long)

part_dros_all_pm <- edgep(temp_out, burnin = burn.in)
rm(temp_out); gc(verbose=FALSE)

write.table(part_dros_all_pm, "par_dros_cont_all_iv_pm.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Order MCMC
temp_out <- orderMCMC(iterations = niter.long,
                          startspace = am_all,
                          scorepar = score_dros,
                          blacklist = blacklist,
                          stepsave = stepsize.long,
                          chainout = TRUE,
                          MAP = FALSE,
                          verbose = TRUE)

# Calculate the posterior probability adjacency matrix.
order_dros_pm <- edgep(temp_out, burnin = burn.in)
rm(temp_out); gc(verbose=FALSE)

write.table(part_dros_pm, "ord_dros_cont_pm.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
