library (baycn)

load('data_geuvadis_pcs.RData')

# Adjacency matrices -----------------------------------------------------------

# adjacency matrix for graphs with five nodes.
am_5 <- matrix(1,
               nrow = 5,
               ncol = 5)

diag(am_5) <- 0

# adjacency matrix for graphs with six nodes.
am_6 <- matrix(1,
               nrow = 6,
               ncol = 6)

diag(am_6) <- 0

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

# Analyze GEUVADIS data with PCs -----------------------------------------------

set.seed(917)

# Run baycn on Q21 with PCs 1 and 3.
bay_Q21_pc <- mhEdge(adjMatrix = am_6,
                     burnIn = 0.2,
                     data = data_Q21_pc,
                     iterations = 50000,
                     nGV = 1,
                     pmr = TRUE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = TRUE,
                     thinTo = 200)

# Run baycn on Q23 with PCs 1 and 8.
bay_Q23_pc <- mhEdge(adjMatrix = am_5,
                     burnIn = 0.2,
                     data = data_Q23_pc,
                     iterations = 50000,
                     nGV = 1,
                     pmr = TRUE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = TRUE,
                     thinTo = 200)

# Run baycn on Q37 with PC 1.
bay_Q37_pc <- mhEdge(adjMatrix = am_5,
                     burnIn = 0.2,
                     data = data_Q37_pc,
                     iterations = 50000,
                     nGV = 1,
                     pmr = TRUE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = TRUE,
                     thinTo = 200)

# Run baycn on Q50 with PCs 1 and 5.
bay_Q50_pc <- mhEdge(adjMatrix = am_5,
                     burnIn = 0.2,
                     data = data_Q50_pc,
                     iterations = 50000,
                     nGV = 1,
                     pmr = TRUE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = TRUE,
                     thinTo = 200)

# Run baycn on Q8 with PCs 1, 2, 6, 7, and 9.
bay_Q8_pc <- mhEdge(adjMatrix = am_9,
                    burnIn = 0.2,
                    data = data_Q8_pc,
                    iterations = 500000,
                    nGV = 1,
                    pmr = TRUE,
                    prior = c(0.05,
                              0.05,
                              0.9),
                    progress = TRUE,
                    thinTo = 1000)

q8_baycn_500k_1k_1 <- bay_Q8_pc@posteriorES

bay_Q8_pc <- mhEdge(adjMatrix = am_9,
                    burnIn = 0.2,
                    data = data_Q8_pc,
                    iterations = 500000,
                    nGV = 1,
                    pmr = TRUE,
                    prior = c(0.05,
                              0.05,
                              0.9),
                    progress = TRUE,
                    thinTo = 1000)
q8_baycn_500k_1k_2 <- bay_Q8_pc@posteriorES

bay_Q8_pc <- mhEdge(adjMatrix = am_9,
                    burnIn = 0.2,
                    data = data_Q8_pc,
                    iterations = 500000,
                    nGV = 1,
                    pmr = TRUE,
                    prior = c(0.05,
                              0.05,
                              0.9),
                    progress = TRUE,
                    thinTo = 1000)
q8_baycn_500k_1k_3 <- bay_Q8_pc@posteriorES

q8_baycn_500k_1k_merged <- q8_baycn_500k_1k_1
q8_baycn_500k_1k_merged$zero <- (q8_baycn_500k_1k_1$zero + q8_baycn_500k_1k_2$zero + q8_baycn_500k_1k_3$zero)/3
q8_baycn_500k_1k_merged$one <- (q8_baycn_500k_1k_1$one + q8_baycn_500k_1k_2$one + q8_baycn_500k_1k_3$one)/3
q8_baycn_500k_1k_merged$two <- (q8_baycn_500k_1k_1$two + q8_baycn_500k_1k_2$two + q8_baycn_500k_1k_3$two)/3

save(bay_Q21_pc,
     bay_Q23_pc,
     bay_Q37_pc,
     bay_Q50_pc,
     bay_Q8_pc,
     file = 'bay_geuvadis_pcs.RData')

write.table(q8_baycn_500k_1k_merged, "q8_baycn_500k_1k_merged_posteriorES.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
