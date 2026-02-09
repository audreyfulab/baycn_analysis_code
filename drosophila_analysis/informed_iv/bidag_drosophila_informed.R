# Load necessary packages ------------------------------------------------------

library (BiDAG)

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

# Generate blacklist to prevent edges from tissue type to TF
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
                          startspace = am_informed,
                          blacklist = blacklist,
                          stepsave = stepsize.long)

part_dros_pm <- edgep(temp_out, burnin = burn.in)
rm(temp_out); gc(verbose=FALSE)

write.table(part_dros_pm, "par_dros_cont_pm.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Order MCMC
temp_out <- orderMCMC(iterations = niter.long,
                          startspace = am_informed,
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
