# data prepration for PCA
rownames(marix_pareto) <- rownames(sampleinfo)
PCA_file <- list(sampleinfo, metainfo, marix_pareto)
names(PCA_file) <- c("sampleinfo", "metainfo", "marix_pareto")

# PCA
dir.create(paste0(pathload,"/PCA Results"))
pca.path <- paste0(pathload,"/PCA Results")
pca_result <- opls(x = PCA_file$marix_pareto, predI = 2)
# plot(pca_result)
pac_modelscore <- pca_result@modelDF
pac_modelscore$PC = rownames(pac_modelscore)
# write.csv(pac_modelscore, paste0(pca.path,"/[1] PCA modelscore.csv"),
#           row.names = F)

# PCA model score
PCA_modelscore <- ggbarplot(pac_modelscore, x = "PC", y = "R2X",
                            fill = "PC",
                            color = "black",
                            palette = "jco",
                            width = 0.5,
                            font.main = c(16, "plain", "black"),            
                            font.x = c(16, "plain", "black"),                  
                            font.y = c(16, "plain", "black"),                
                            font.legend = c(16, "plain", "black"),
                            label = T,
                            font.label = 16, repel = TRUE,
                            xlab = "Components",
                            ylab = "R2X Value",
                            legend = "top",
                            legend.title = "Components") +
  theme(axis.text = element_text(size = 16))
ggsave("[1] PCA modelscore.png", 
       egg::set_panel_size(PCA_modelscore, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = pca.path)

ggsave("[1] PCA modelscore.pdf", 
       egg::set_panel_size(PCA_modelscore, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = pca.path)

# PCA score plot
scoreMN <- pca_result@scoreMN
scoreMN <- cbind(scoreMN, PCA_file$sampleinfo)
scoreMN$samples <- rownames(scoreMN)
# write.csv(scoreMN,paste0(pca.path,"/[2] PCA score.csv"))


PCA_score <- ggscatter(scoreMN, x = "p1", y = "p2",
                       color = "black", 
                       fill = "Treatment",
                       shape = 21,
                       ellipse = T, 
                       size = 5,
                       ellipse.level = 0.75,
                       mean.point = F,
                       alpha = 0.8,
                       font.label = 14, repel = TRUE,
                       xlab = paste0("PC1 (",pca_result@modelDF$R2X[1]*100,"%)"),
                       ylab = paste0("PC2 (",pca_result@modelDF$R2X[2]*100,"%)"),
                       font.main = c(16, "plain", "black"),            
                       font.x = c(16, "plain", "black"),                  
                       font.y = c(16, "plain", "black"),                
                       legend = "top",  
                       legend.title = "Treatment",  
                       font.legend = c(16, "plain", "black"),               
                       rotate = F,                                         
                       ticks = T,   
                       label = "samples",
                       tickslab = T,  
                       palette = "jco") +
  theme(axis.text = element_text(size = 16))

ggsave("[2] PCA score.png",
       egg::set_panel_size(PCA_score, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = pca.path)
ggsave("[2] PCA score.pdf",
       egg::set_panel_size(PCA_score, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = pca.path)


# PCA loading plot
loadingMN <- pca_result@loadingMN
loadingMN <- as.data.frame(loadingMN)
loadingMN$compound <- as.vector(rownames(loadingMN))
colnames(loadingMN)[1:2] = c("PC1","PC2")

loadingMN_sort <- loadingMN %>% 
  mutate(distance = sqrt(PC1^2+PC2^2)) %>%
  arrange(desc(distance))

loadingMN_label <- loadingMN_sort %>%
  mutate(label = ifelse(distance > loadingMN_sort$distance[20], compound, "")) %>%
  mutate(alphavalue = ifelse(distance > loadingMN_sort$distance[20], 1, 0.8))

PCA_loading <- ggscatter(loadingMN_label, x = "PC1", y = "PC2",
                         fill = "distance", 
                         shape = 21,
                         size =  3,
                         alpha = "alphavalue",
                         legend = " ",
                         font.main = c(16, "plain", "black"),            
                         font.x = c(16, "plain", "black"),                  
                         font.y = c(16, "plain", "black"),                
                         font.legend = c(16, "plain", "black"),
                         xlab = "Loading 1",
                         ylab = "Loading 2") +
  scale_color_gradient(low = "grey", high = "#0077C2",na.value = NA) +
  geom_text_repel(data = loadingMN_label, size = 4, aes(label = label))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text = element_text(size = 16)) 
ggsave("[3] PCA loading.png", 
       egg::set_panel_size(PCA_loading, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = pca.path)
ggsave("[3] PCA loading.pdf", 
       egg::set_panel_size(PCA_loading, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = pca.path)


wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
addWorksheet(wb, "Sheet 3", tabColour = "#008000")

writeData(wb, sheet = 1, pac_modelscore)
writeData(wb, sheet = 2, scoreMN)
writeData(wb, sheet = 3, loadingMN_label)

names(wb)[[1]] <- "PCA modelscore"
names(wb)[[2]] <- "PCA score"
names(wb)[[3]] <- "PCA loading"

saveWorkbook(wb, paste0(pca.path,"/pca data.xlsx"), overwrite = TRUE)

