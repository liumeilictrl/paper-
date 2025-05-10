rm(list = ls())
library(magrittr)
library(tidygraph)
library(ggraph)
library(psych)
library(tibble)
dat <- fread("output/LUSC_TPM.csv", data.table = F)
dat <- dat %>% column_to_rownames(var ='V1')
genes <- fread("input/Pyroptosis related genes.csv", data.table = F)
genes <- genes$genes
dat <- dat[genes,] %>% t()

res_cor <- corr.test(x = dat, y = dat, method = "spearman", adjust = "fdr")
res_cor_r <- reshape::melt(res_cor$r)
res_cor_p <- reshape::melt(res_cor$p.adj)
res_cor_all <- cbind(res_cor_r, res_cor_p[,3])
colnames(res_cor_all) <- c("gene1", "gene2", "cor", "p.adj")
write.csv(res_cor_all, "output/res_cor.csv")

edges <- res_cor_r
colnames(edges) <- c("from", "to", "cor")
nodes <- data.frame(id = genes)
net.tidy <- tbl_graph(edges = edges, nodes = nodes, directed = TRUE)

set.seed(123)
ggraph(net.tidy, layout = "graphopt") + 
  geom_edge_link(aes(width = abs(cor), alpha = abs(cor), color = cor)) +
  scale_edge_width(range = c(0.2, 2)) +
  scale_edge_color_gradient2(low = "blue", high = "red", mid = "white") +
  geom_node_point(colour = "white" , size = 6) +
  geom_node_text(aes(label = id), repel = TRUE, check_overlap = T, size = 5) + 
  theme_graph()

ggsave("output/network_cor.pdf", width = 8, height = 8, device = cairo_pdf)
