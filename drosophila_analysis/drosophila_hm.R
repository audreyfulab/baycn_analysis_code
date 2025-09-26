# Load libraries ---------------------------------------------------------------

library (reshape2)
library (ggplot2)
library (baycn)

# Load drosophila data ---------------------------------------------------------

# Read in the drosophila data.
#drosophila_raw <- read.delim('/Users/Evatar/Downloads/drosophila.txt',
#                             comment.char = '#')

# Remove the rows with genetic information.
#drosophila_continuous <- drosophila_raw[, -c(1:5)]

# Load the drosphila data from baycn.
data("drosophila")
length (drosophila)
# 2
summary (drosophila)
#Length Class  Mode   
#continuous 6510   -none- numeric
#discrete   6510   -none- numeric
dim (drosophila$discrete)
# 310  21
# the first six columns are tissues

# Correlation heatmap: discrete ------------------------------------------------

# Calculate the correlation of the discrete drosophila data.
cor_dros_dis <- round(cor(drosophila$discrete), 2)

# Melt the correlation matrix.
melted_cor_dros_dis <- melt(cor_dros_dis)

# plot the correlation heatmap.
p_dros_dis <- ggplot(data = melted_cor_dros_dis,
                     aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  ggtitle('') +
  # geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  scale_fill_gradient2(low = 'darkblue', high = 'darkorange', mid = 'white',
                       midpoint = 0, limit = c(-1, 1), space = 'Lab',
                       name = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        legend.position = 'right')

pdf ("Plot_heatmap_drosophila_dis_cor.pdf")
p_dros_dis
dev.off()
# ggsave(filename = 'drosophila_hm.pdf',
#        plot = p_dros_dis,
#        height = 7 ,
#        width = 7,
#        path = '/Users/Evatar/Desktop/')

# Correlation heatmap: continuous ----------------------------------------------

# Calculate the correlation of the discrete drosophila data.
cor_dros_con <- round(cor(drosophila$continuous), 2)

# Melt the correlation matrix.
melted_cor_dros_con <- melt(cor_dros_con)

# plot the correlation heatmap.
p_dros_con <- ggplot(data = melted_cor_dros_con,
                     aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  ggtitle('') +
  # geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  scale_fill_gradient2(low = 'darkblue', high = 'darkorange', mid = 'white',
                       midpoint = 0, limit = c(-1, 1), space = 'Lab',
                       name = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        legend.position = 'right')

pdf ("Plot_heatmap_drosophila_con_cor.pdf")
p_dros_con
dev.off()

# ggsave(filename = 'drosophila_hm_con.pdf',
#        plot = p_dros_con,
#        height = 7 ,
#        width = 7,
#        path = '/Users/Evatar/Desktop/')
