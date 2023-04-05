dir.create(paste0(pathload,"/Cluster"))
cluster.path <- paste0(pathload,"/Cluster")
# 
# BiocManager::install("ggtree")
# library(ggtree)


data_gapfilled$Treatment %>% levels() %>% length() -> k

USArrests <- data_gapfilled %>%
  select(-Treatment) %>%
  column_to_rownames( var = "Sample") 

dd <- dist(scale(USArrests), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

print(factoextra::fviz_dend(hc, cex = 1, k_colors = "jco",
                            main = " ",k = k,
                            
                            xlab = " ", ylab = "Distance", sub = "")) 



pdf(paste0(cluster.path,"/[1] Clustering.pdf"), height = 4, width = 5)
print(factoextra::fviz_dend(hc, cex = 0.8, k_colors = "jco",
                            main = " ",k = k,
                            xlab = " ", ylab = "Distance", sub = ""))
dev.off()


png(paste0(cluster.path,"/[1] Clustering.png"),height = 1000, width = 1300,
    units = "px",res = 256)
print(factoextra::fviz_dend(hc, cex = 0.8, k_colors = "jco",
                            main = " ",k = k,
                            xlab = " ", ylab = "Distance", sub = ""))
dev.off()


pdf(paste0(cluster.path,"/[2] Clustering.pdf"),height = 4, width = 5)
print(factoextra::fviz_dend(hc, cex = 0.6, k = k,
                            k_colors = "jco", type = "circular"))
dev.off()


png(paste0(cluster.path,"/[2] Clustering.png"),height = 1000, width = 1000,
    units = "px",res = 256)
print(factoextra::fviz_dend(hc, cex = 0.6, k = k,
                            k_colors = "jco", type = "circular"))
dev.off()


