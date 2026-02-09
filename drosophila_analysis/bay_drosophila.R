# Load the drosophila data -----------------------------------------------------

setwd("~/Documents/Programs/baycn_analysis_code/drosophila_analysis")
load('data_drosophila_tissue.RData')

# Adjacency matrices -----------------------------------------------------------
am_dros <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  byrow = TRUE,
                  nrow = 16)

# Run baycn on each tissue type separately -------------------------------------

# baycn: Meso --------------------------

set.seed(22)

baycn_meso <- mhEdge(data = data_meso,
                     adjMatrix = am_dros,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     nCPh = 1,
                     nGV = 0,
                     pmr = FALSE,
                     burnIn = 0.2,
                     iterations = 30000,
                     thinTo = 500,
                     progress = TRUE)

baycn_meso_dis <- mhEdge(data = data_dis_meso,
                         adjMatrix = am_dros,
                         prior = c(0.05,
                                   0.05,
                                   0.9),
                         nCPh = 1,
                         nGV = 0,
                         pmr = FALSE,
                         burnIn = 0.2,
                         iterations = 30000,
                         thinTo = 500,
                         progress = TRUE)

# baycn: VM ----------------------------

set.seed(919)

baycn_vm <- mhEdge(data = data_vm,
                   adjMatrix = am_dros,
                   prior = c(0.05,
                             0.05,
                             0.9),
                   nCPh = 1,
                   nGV = 0,
                   pmr = FALSE,
                   burnIn = 0.2,
                   iterations = 30000,
                   thinTo = 500,
                   progress = TRUE)

baycn_vm_dis <- mhEdge(data = data_dis_vm,
                       adjMatrix = am_dros,
                       prior = c(0.05,
                                 0.05,
                                 0.9),
                       nCPh = 1,
                       nGV = 0,
                       pmr = FALSE,
                       burnIn = 0.2,
                       iterations = 30000,
                       thinTo = 500,
                       progress = TRUE)

# baycn: SM ----------------------------

set.seed(228)

baycn_sm <- mhEdge(data = data_sm,
                   adjMatrix = am_dros,
                   prior = c(0.05,
                             0.05,
                             0.9),
                   nCPh = 1,
                   nGV = 0,
                   pmr = FALSE,
                   burnIn = 0.2,
                   iterations = 30000,
                   thinTo = 500,
                   progress = TRUE)

baycn_sm_dis <- mhEdge(data = data_dis_sm,
                       adjMatrix = am_dros,
                       prior = c(0.05,
                                 0.05,
                                 0.9),
                       nCPh = 1,
                       nGV = 0,
                       pmr = FALSE,
                       burnIn = 0.2,
                       iterations = 30000,
                       thinTo = 500,
                       progress = TRUE)

# baycn: CM ----------------------------

set.seed(669)

baycn_cm <- mhEdge(data = data_cm,
                   adjMatrix = am_dros,
                   prior = c(0.05,
                             0.05,
                             0.9),
                   nCPh = 1,
                   nGV = 0,
                   pmr = FALSE,
                   burnIn = 0.2,
                   iterations = 30000,
                   thinTo = 500,
                   progress = TRUE)

baycn_cm_dis <- mhEdge(data = data_dis_cm,
                       adjMatrix = am_dros,
                       prior = c(0.05,
                                 0.05,
                                 0.9),
                       nCPh = 1,
                       nGV = 0,
                       pmr = FALSE,
                       burnIn = 0.2,
                       iterations = 30000,
                       thinTo = 500,
                       progress = TRUE)

# baycn: Meso_SM -----------------------

set.seed(231)

baycn_meso_sm <- mhEdge(data = data_meso_sm,
                        adjMatrix = am_dros,
                        prior = c(0.05,
                                  0.05,
                                  0.9),
                        nCPh = 1,
                        nGV = 0,
                        pmr = FALSE,
                        burnIn = 0.2,
                        iterations = 30000,
                        thinTo = 500,
                        progress = TRUE)

baycn_meso_sm_dis <- mhEdge(data = data_dis_meso_sm,
                            adjMatrix = am_dros,
                            prior = c(0.05,
                                      0.05,
                                      0.9),
                            nCPh = 1,
                            nGV = 0,
                            pmr = FALSE,
                            burnIn = 0.2,
                            iterations = 30000,
                            thinTo = 500,
                            progress = TRUE)

# baycn: VM_SM -------------------------

set.seed(396)

baycn_vm_sm <- mhEdge(data = data_vm_sm,
                      adjMatrix = am_dros,
                      prior = c(0.05,
                                0.05,
                                0.9),
                      nCPh = 1,
                      nGV = 0,
                      pmr = FALSE,
                      burnIn = 0.2,
                      iterations = 30000,
                      thinTo = 500,
                      progress = TRUE)

baycn_vm_sm_dis <- mhEdge(data = data_dis_vm_sm,
                          adjMatrix = am_dros,
                          prior = c(0.05,
                                    0.05,
                                    0.9),
                          nCPh = 1,
                          nGV = 0,
                          pmr = FALSE,
                          burnIn = 0.2,
                          iterations = 30000,
                          thinTo = 500,
                          progress = TRUE)

# Plot the inferred graph for each tissue --------------------------------------

set.seed(8)
igraph_options(plot.layout=layout_nicely)
plot(baycn_meso, vertex.shape = 'none', edge.arrow.size = 0.25)
plot(baycn_vm, vertex.shape = 'none', edge.arrow.size = 0.25)
plot(baycn_sm, vertex.shape = 'none', edge.arrow.size = 0.25)
plot(baycn_cm, vertex.shape = 'none', edge.arrow.size = 0.25)
plot(baycn_meso_sm, vertex.shape = 'none', edge.arrow.size = 0.25)
plot(baycn_vm_sm, vertex.shape = 'none', edge.arrow.size = 0.25)

# Save the output from baycn for each tissue -----------------------------------

# save(baycn_meso,
#      baycn_vm,
#      baycn_sm,
#      baycn_cm,
#      baycn_meso_sm,
#      baycn_vm_sm,
#      baycn_meso_dis,
#      baycn_vm_dis,
#      baycn_sm_dis,
#      baycn_cm_dis,
#      baycn_meso_sm_dis,
#      baycn_vm_sm_dis,
#      file = 'bay_drosophila.RData')

# Run lasso on each tissue -----------------------------------------------------

set.seed(608)

lasso_meso <- cv.glmnet(x = data_meso[, 1:15],
                        y = data_meso[, 16],
                        family = 'binomial',
                        alpha = 1)

coef(lasso_meso,
     s = lasso_meso$lambda.min)

lasso_vm <- cv.glmnet(x = data_vm[, 1:15],
                      y = data_vm[, 16],
                      family = 'binomial',
                      alpha = 1)

coef(lasso_vm,
     s = lasso_vm$lambda.min)

lasso_sm <- cv.glmnet(x = data_sm[, 1:15],
                      y = data_sm[, 16],
                      family = 'binomial',
                      alpha = 1)

coef(lasso_sm,
     s = lasso_sm$lambda.min)

lasso_cm <- cv.glmnet(x = data_cm[, 1:15],
                      y = data_cm[, 16],
                      family = 'binomial',
                      alpha = 1)

coef(lasso_cm,
     s = lasso_cm$lambda.min)

lasso_meso_sm <- cv.glmnet(x = data_meso_sm[, 1:15],
                           y = data_meso_sm[, 16],
                           family = 'binomial',
                           alpha = 1)

coef(lasso_meso_sm,
     s = lasso_meso_sm$lambda.1se)

lasso_vm_sm <- cv.glmnet(x = data_vm_sm[, 1:15],
                         y = data_vm_sm[, 16],
                         family = 'binomial',
                         alpha = 1)

coef(lasso_vm_sm,
     s = lasso_vm_sm$lambda.min)
