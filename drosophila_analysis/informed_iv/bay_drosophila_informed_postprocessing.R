# process baycn output from long runs

# runs of 5m iterations
load("dros_baycn_5e+06_1.RData")
baycn_dros_posteriorES_1 <- baycn_dros@posteriorES

load("dros_baycn_5e+06_2.RData")
baycn_dros_posteriorES_2 <- baycn_dros@posteriorES

load("dros_baycn_5e+06_3.RData")
baycn_dros_posteriorES_3 <- baycn_dros@posteriorES

dros_baycn_5m_1k_merged <- baycn_dros_posteriorES_1
dros_baycn_5m_1k_merged$zero <- (baycn_dros_posteriorES_1$zero + baycn_dros_posteriorES_2$zero + baycn_dros_posteriorES_3$zero)/3
dros_baycn_5m_1k_merged$one <- (baycn_dros_posteriorES_1$one + baycn_dros_posteriorES_2$one + baycn_dros_posteriorES_3$one)/3
dros_baycn_5m_1k_merged$two <- (baycn_dros_posteriorES_1$two + baycn_dros_posteriorES_2$two + baycn_dros_posteriorES_3$two)/3

write.table(dros_baycn_5m_1k_merged, "dros_baycn_informed_iv_5m_1k_merged_posteriorES.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

dros_baycn_5m_1k_merged[grep("-Meso", dros_baycn_5m_1k_merged$nodes),]

# runs of 50m iterations
load("dros_baycn_5e+07_1.RData")
baycn_dros_posteriorES_1 <- baycn_dros@posteriorES

load("dros_baycn_5e+07_2.RData")
baycn_dros_posteriorES_2 <- baycn_dros@posteriorES

load("dros_baycn_5e+07_3.RData")
baycn_dros_posteriorES_3 <- baycn_dros@posteriorES

dros_baycn_50m_1k_merged <- baycn_dros_posteriorES_1
dros_baycn_50m_1k_merged$zero <- (baycn_dros_posteriorES_1$zero + baycn_dros_posteriorES_2$zero + baycn_dros_posteriorES_3$zero)/3
dros_baycn_50m_1k_merged$one <- (baycn_dros_posteriorES_1$one + baycn_dros_posteriorES_2$one + baycn_dros_posteriorES_3$one)/3
dros_baycn_50m_1k_merged$two <- (baycn_dros_posteriorES_1$two + baycn_dros_posteriorES_2$two + baycn_dros_posteriorES_3$two)/3

write.table(dros_baycn_50m_1k_merged, "dros_baycn_50m_1k_merged_posteriorES.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

dros_baycn_50m_1k_merged[grep("-Meso", dros_baycn_50m_1k_merged$nodes),]
