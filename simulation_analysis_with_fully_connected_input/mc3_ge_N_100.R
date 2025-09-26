#install.packages("remotes")
#remotes::install_github("rjbgoudie/parental")
#remotes::install_github("rjbgoudie/structmcmc")

# structmcmc calls .Internal(crossprod(...))
# Since R ≥4.3 made crossprod() a primitive instead of an .Internal, any such call now fails
# 
# Need the following patch:
#
# patch_crossprod_calls <- function(pkg){
#   ns <- asNamespace(pkg)
#   fixed <- character()
#   for (fname in ls(ns, all.names = TRUE)){
#     f <- get(fname, envir = ns, inherits = FALSE)
#     if (!is.function(f)) next
#     bod <- paste(deparse(body(f)), collapse = "\n")
#     if (!grepl("\\.Internal\\(\\s*crossprod\\s*\\(", bod)) next
#     
#     bod2 <- bod
#     bod2 <- gsub("\\.Internal\\(\\s*crossprod\\s*\\(([^,\\)]+),\\s*NULL\\s*\\)\\)", "base::crossprod(\\1)", bod2)
#     bod2 <- gsub("\\.Internal\\(\\s*crossprod\\s*\\(([^\\)]+)\\)\\)", "base::crossprod(\\1)", bod2)
#     
#     # apply patch
#     f2 <- f
#     body(f2) <- parse(text = bod2)[[1]]
#     unlockBinding(fname, ns); assign(fname, f2, envir = ns); lockBinding(fname, ns)
#     fixed <- c(fixed, paste(pkg, fname, sep="::"))
#   }
#   fixed
# }
# 
# fixed <- c(patch_crossprod_calls("parental"),
#            patch_crossprod_calls("structmcmc"))
# fixed


library (parental)
library (structmcmc)

load('./data_ge_N_100.RData')

# Constrained matricies were all removed. Default behavior is Fully Connected

# ------------------------------------

# Initialize each list for the gibbs output to full length.
mc3_g2_100_02_fc <- vector(mode = 'list',
                           length = M)

mc3_g2_100_05_fc <- vector(mode = 'list',
                           length = M)

mc3_g2_100_1_fc <- vector(mode = 'list',
                          length = M)

mc3_g2_100_02 <- vector(mode = 'list',
                        length = M)

mc3_g2_100_05 <- vector(mode = 'list',
                        length = M)

mc3_g2_100_1 <- vector(mode = 'list',
                       length = M)

mc3_nc11_100_02 <- vector(mode = 'list',
                          length = M)

mc3_nc11_100_05 <- vector(mode = 'list',
                          length = M)

mc3_nc11_100_1 <- vector(mode = 'list',
                         length = M)

mc3_pc_100_02 <- vector(mode = 'list',
                        length = M)

mc3_pc_100_05 <- vector(mode = 'list',
                        length = M)

mc3_pc_100_1 <- vector(mode = 'list',
                       length = M)

# Initialize the probability matrix lists to full length.
mc3_g2_100_02_fc_pm <- vector(mode = 'list',
                              length = M)

mc3_g2_100_05_fc_pm <- vector(mode = 'list',
                              length = M)

mc3_g2_100_1_fc_pm <- vector(mode = 'list',
                             length = M)

mc3_g2_100_02_pm <- vector(mode = 'list',
                           length = M)

mc3_g2_100_05_pm <- vector(mode = 'list',
                           length = M)

mc3_g2_100_1_pm <- vector(mode = 'list',
                          length = M)

mc3_nc11_100_02_pm <- vector(mode = 'list',
                             length = M)

mc3_nc11_100_05_pm <- vector(mode = 'list',
                             length = M)

mc3_nc11_100_1_pm <- vector(mode = 'list',
                            length = M)

mc3_pc_100_02_pm <- vector(mode = 'list',
                           length = M)

mc3_pc_100_05_pm <- vector(mode = 'list',
                           length = M)

mc3_pc_100_1_pm <- vector(mode = 'list',
                          length = M)

# Initialize the time lists to full length.
mc3_g2_100_02_fc_time <- vector(mode = 'list',
                                length = M)

mc3_g2_100_05_fc_time <- vector(mode = 'list',
                                length = M)

mc3_g2_100_1_fc_time <- vector(mode = 'list',
                               length = M)

mc3_g2_100_02_time <- vector(mode = 'list',
                             length = M)

mc3_g2_100_05_time <- vector(mode = 'list',
                             length = M)

mc3_g2_100_1_time <- vector(mode = 'list',
                            length = M)

mc3_nc11_100_02_time <- vector(mode = 'list',
                               length = M)

mc3_nc11_100_05_time <- vector(mode = 'list',
                               length = M)

mc3_nc11_100_1_time <- vector(mode = 'list',
                              length = M)

mc3_pc_100_02_time <- vector(mode = 'list',
                             length = M)

mc3_pc_100_05_time <- vector(mode = 'list',
                             length = M)

mc3_pc_100_1_time <- vector(mode = 'list',
                            length = M)

set.seed(101)

# Loop through each combination of signal strength and sample size for all
# topologies.
for (e in 1:M) {

  #################################################
  # G2 FULLY CONNECTED
  #################################################

  starttime <- Sys.time()
  mc3_g2_100_02_fc[[e]] <- posterior(data = data_g2_100_02[[e]],
                                     method = 'mc3',
                                     logScoreFUN = logScoreNormalFUN(),
                                     nSamples = 30000,
                                     nBurnin = 6000,
                                     verbose = FALSE)
  endtime <- Sys.time()
  mc3_g2_100_02_fc_time[[e]] <- endtime - starttime
  mc3_g2_100_02_fc_pm[[e]] <- ep(mc3_g2_100_02_fc[[e]])

  starttime <- Sys.time()
  mc3_g2_100_05_fc[[e]] <- posterior(data = data_g2_100_05[[e]],
                                     method = 'mc3',
                                     logScoreFUN = logScoreNormalFUN(),
                                     nSamples = 30000,
                                     nBurnin = 6000,
                                     verbose = FALSE)
  endtime <- Sys.time()
  mc3_g2_100_05_fc_time[[e]] <- endtime - starttime
  mc3_g2_100_05_fc_pm[[e]] <- ep(mc3_g2_100_05_fc[[e]])

  starttime <- Sys.time()
  mc3_g2_100_1_fc[[e]] <- posterior(data = data_g2_100_1[[e]],
                                    method = 'mc3',
                                    logScoreFUN = logScoreNormalFUN(),
                                    nSamples = 30000,
                                    nBurnin = 6000,
                                    verbose = FALSE)
  endtime <- Sys.time()
  mc3_g2_100_1_fc_time[[e]] <- endtime - starttime
  mc3_g2_100_1_fc_pm[[e]] <- ep(mc3_g2_100_1_fc[[e]])

  print(Sys.time())

  #################################################
  # NC11
  #################################################

  starttime <- Sys.time()
  mc3_nc11_100_02[[e]] <- posterior(data = data_nc11_100_02[[e]],
                                    method = 'mc3',
                                    logScoreFUN = logScoreNormalFUN(),
                                    nSamples = 50000,
                                    nBurnin = 10000,
                                    verbose = FALSE)
  endtime <- Sys.time()
  mc3_nc11_100_02_time[[e]] <- endtime - starttime
  mc3_nc11_100_02_pm[[e]] <- ep(mc3_nc11_100_02[[e]])

  starttime <- Sys.time()
  mc3_nc11_100_05[[e]] <- posterior(data = data_nc11_100_05[[e]],
                                    method = 'mc3',
                                    logScoreFUN = logScoreNormalFUN(),
                                    nSamples = 50000,
                                    nBurnin = 10000,
                                    verbose = FALSE)
  endtime <- Sys.time()
  mc3_nc11_100_05_time[[e]] <- endtime - starttime
  mc3_nc11_100_05_pm[[e]] <- ep(mc3_nc11_100_05[[e]])

  starttime <- Sys.time()
  mc3_nc11_100_1[[e]] <- posterior(data = data_nc11_100_1[[e]],
                                   method = 'mc3',
                                   logScoreFUN = logScoreNormalFUN(),
                                   nSamples = 50000,
                                   nBurnin = 10000,
                                   verbose = FALSE)
  endtime <- Sys.time()
  mc3_nc11_100_1_time[[e]] <- endtime - starttime
  mc3_nc11_100_1_pm[[e]] <- ep(mc3_nc11_100_1[[e]])

  print('nc11')
  print(Sys.time())

  #################################################
  # PC
  #################################################

  starttime <- Sys.time()
  mc3_pc_100_02[[e]] <- posterior(data = data_pc_100_02[[e]],
                                  method = 'mc3',
                                  logScoreFUN = logScoreNormalFUN(),
                                  nSamples = 50000,
                                  nBurnin = 10000,
                                  verbose = FALSE)
  endtime <- Sys.time()
  mc3_pc_100_02_time[[e]] <- endtime - starttime
  mc3_pc_100_02_pm[[e]] <- ep(mc3_pc_100_02[[e]])

  starttime <- Sys.time()
  mc3_pc_100_05[[e]] <- posterior(data = data_pc_100_05[[e]],
                                  method = 'mc3',
                                  logScoreFUN = logScoreNormalFUN(),
                                  nSamples = 50000,
                                  nBurnin = 10000,
                                  verbose = FALSE)
  endtime <- Sys.time()
  mc3_pc_100_05_time[[e]] <- endtime - starttime
  mc3_pc_100_05_pm[[e]] <- ep(mc3_pc_100_05[[e]])

  starttime <- Sys.time()
  mc3_pc_100_1[[e]] <- posterior(data = data_pc_100_1[[e]],
                                 method = 'mc3',
                                 logScoreFUN = logScoreNormalFUN(),
                                 nSamples = 50000,
                                 nBurnin = 10000,
                                 verbose = FALSE)
  endtime <- Sys.time()
  mc3_pc_100_1_time[[e]] <- endtime - starttime
  mc3_pc_100_1_pm[[e]] <- ep(mc3_pc_100_1[[e]])

  print('pc')
  print(e)
  print(Sys.time())

  save(M,
       mc3_g2_100_02_fc,
       mc3_g2_100_05_fc,
       mc3_g2_100_1_fc,
       mc3_g2_100_02,
       mc3_g2_100_05,
       mc3_g2_100_1,
       mc3_nc11_100_02,
       mc3_nc11_100_05,
       mc3_nc11_100_1,
       mc3_pc_100_02,
       mc3_pc_100_05,
       mc3_pc_100_1,
       file = './mc3_ge_N_100.RData')

  save(M,
       mc3_g2_100_02_fc_pm,
       mc3_g2_100_05_fc_pm,
       mc3_g2_100_1_fc_pm,
       mc3_g2_100_02_pm,
       mc3_g2_100_05_pm,
       mc3_g2_100_1_pm,
       mc3_nc11_100_02_pm,
       mc3_nc11_100_05_pm,
       mc3_nc11_100_1_pm,
       mc3_pc_100_02_pm,
       mc3_pc_100_05_pm,
       mc3_pc_100_1_pm,
       file = './mc3_ge_N_100_pm.RData')

  save(M,
       mc3_g2_100_02_fc_time,
       mc3_g2_100_05_fc_time,
       mc3_g2_100_1_fc_time,
       mc3_g2_100_02_time,
       mc3_g2_100_05_time,
       mc3_g2_100_1_time,
       mc3_nc11_100_02_time,
       mc3_nc11_100_05_time,
       mc3_nc11_100_1_time,
       mc3_pc_100_02_time,
       mc3_pc_100_05_time,
       mc3_pc_100_1_time,
       file = './mc3_ge_N_100_time.RData')

}
