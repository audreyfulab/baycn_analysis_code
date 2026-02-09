library(baycn)

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
                     iterations = 5e4,
                     thinTo = 4e4,
                     progress = TRUE)

svg("Trace_loglik_graph_dros_informed_50k_1.svg")
tracePlot(baycn_dros)
dev.off()

baycn_dros <- mhEdge(data = drosophila$continuous,
                     adjMatrix = am_informed,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     nCPh = 0,
                     nGV = 6,
                     pmr = TRUE,
                     burnIn = 0.2,
                     iterations = 5e4,
                     thinTo = 4e4,
                     progress = TRUE)

svg("Trace_loglik_graph_dros_informed_50k_2.svg")
tracePlot(baycn_dros)
dev.off()

baycn_dros <- mhEdge(data = drosophila$continuous,
                     adjMatrix = am_informed,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     nCPh = 0,
                     nGV = 6,
                     pmr = TRUE,
                     burnIn = 0.2,
                     iterations = 5e4,
                     thinTo = 4e4,
                     progress = TRUE)

svg("Trace_loglik_graph_dros_informed_50k_3.svg")
tracePlot(baycn_dros)
dev.off()
