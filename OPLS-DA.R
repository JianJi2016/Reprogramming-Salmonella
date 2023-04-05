# OPLS-DA
dir.create(paste0(pathload,"/OPLS-DA Results"))
oplsda.path <- paste0(pathload,"/OPLS-DA Results")

oplsda_result <- opls(x = PCA_file$marix_pareto, 
                      y = PCA_file$sampleinfo[, 'Treatment'],
                      predI = 1, orthoI = 1)
oplsda_modelscore <- oplsda_result@modelDF

oplsda_plot <- oplsda_modelscore[1:2,] %>% 
  tibble::rownames_to_column(var = "group") %>%
  select(group,`R2Y(cum)`,`Q2(cum)`) %>%
  reshape2::melt(.,1)

OPLSDA_modelscore <- ggbarplot(oplsda_plot, x= "group", y="value",
                               fill = "variable",palette = "jco",
                               color = "black",
                               font.main = c(16, "plain", "black"),
                               font.x = c(16, "plain", "black"),
                               font.y = c(16, "plain", "black"),
                               font.legend = c(16, "plain", "black"),
                               label = T,
                               xlab = "Components",
                               ylab = "Score",
                               # legend = "right",
                               legend.title = " ",
                               position = position_dodge(0.7))+
  theme(axis.text = element_text(size = 16))
ggsave("[1] OPLS-DA Model score.png",
       egg::set_panel_size(OPLSDA_modelscore, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = oplsda.path)
ggsave("[1] OPLS-DA Model score.pdf",
       egg::set_panel_size(OPLSDA_modelscore, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = oplsda.path)

# PLSDA score plot
oplsda_scoreMN <- cbind(oplsda_result@scoreMN, oplsda_result@orthoScoreMN)

oplsda_scoreMN <- cbind(oplsda_scoreMN, PCA_file$sampleinfo)
oplsda_scoreMN$samples <- rownames(oplsda_scoreMN)
# write.csv(oplsda_scoreMN,paste0(oplsda.path,"/[2] OPLS-DA score.csv"))

OPLSDA_score <- ggscatter(oplsda_scoreMN, x = "p1", y = "o1",
                          color = "black", 
                          fill = "Treatment",
                          shape = 21,
                          ellipse = T, 
                          size = 5,
                          ellipse.level = 0.75,
                          mean.point = F,
                          alpha = 0.8,
                          font.label = 14, repel = TRUE,
                          xlab = paste0("PC1 (",oplsda_result@modelDF$R2X[1]*100,"%)"),
                          ylab = paste0("PC2 (",oplsda_result@modelDF$R2X[2]*100,"%)"),
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
ggsave("[2] OPLS-DA score.png",
       egg::set_panel_size(OPLSDA_score, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = oplsda.path)
ggsave("[2] OPLS-DA score.pdf",
       egg::set_panel_size(OPLSDA_score, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = oplsda.path)

# OPLSDA loading plot
oplsda_loadingMN <- cbind(oplsda_result@loadingMN,oplsda_result@orthoLoadingMN)
oplsda_loadingMN <- as.data.frame(oplsda_loadingMN)
oplsda_loadingMN$compound <- as.vector(rownames(oplsda_loadingMN))
colnames(oplsda_loadingMN)[1:2] = c("PC1","PC2")

oplsda_loadingMN_sort <- oplsda_loadingMN %>% 
  mutate(distance = PC1^2+PC2^2) %>%
  arrange(desc(distance))

oplsda_loadingMN_label <- oplsda_loadingMN_sort %>%
  mutate(label = ifelse(distance > oplsda_loadingMN_sort$distance[20], compound, "")) %>%
  mutate(alphavalue = ifelse(distance > oplsda_loadingMN_sort$distance[20], 1, 0.8))
# write.csv(oplsda_loadingMN_label, paste0(oplsda.path,"/[3] OPLSDA loading.csv"))

OPLSDA_loading <- ggscatter(oplsda_loadingMN_label, x = "PC1", y = "PC2",
                            fill = "distance", 
                            shape = 21,
                            size =  4,
                            alpha = "alphavalue",
                            legend = " ",
                            font.main = c(16, "plain", "black"),            
                            font.x = c(16, "plain", "black"),                  
                            font.y = c(16, "plain", "black"),                
                            font.legend = c(16, "plain", "black"),
                            xlab = "Loading 1",
                            ylab = "Loading 2") +
  scale_color_gradient(low = "grey", high = "#0077C2") +
  geom_text_repel(data = oplsda_loadingMN_label, size = 4,aes(label = label))+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text = element_text(size = 16))
ggsave("[3] OPLSDA loading.png", 
       egg::set_panel_size(OPLSDA_loading, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = oplsda.path)
ggsave("[3] OPLSDA loading.pdf", 
       egg::set_panel_size(OPLSDA_loading, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = oplsda.path)


# VIP
getVipVn(oplsda_result) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "name") ->medi2
colnames(medi2)[2] = "VIP"
oplsda_VIP <- medi2

# write.csv(oplsda_VIP,paste0(oplsda.path,"/[4] OPLS-DA VIP.csv"))

oplsda_VIPplot <- oplsda_VIP %>% arrange(desc(VIP)) %>%
  slice(1:10)

OLSDA_VIP <- ggbarplot(oplsda_VIPplot, x = "name", y = "VIP",
                       xlab = "",
                       ylab = "VIP Value",
                       fill = "VIP",
                       color = "black",
                       legend = "",
                       font.main = c(16, "plain", "black"),            
                       font.x = c(16, "plain", "black"),                  
                       font.y = c(16, "plain", "black"),                
                       font.legend = c(16, "plain", "black"),
                       orientation = "horiz",
                       order = rev(oplsda_VIPplot$name)) +
  scale_fill_gradient(low = "grey", high = "#0077C2") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text(size = 16))

ggsave("[4] OPLS-DA VIP.png",
       egg::set_panel_size(OLSDA_VIP, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height =5, width = 10, path = oplsda.path)
ggsave("[4] OPLS-DA VIP.pdf",
       egg::set_panel_size(OLSDA_VIP, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height =5, width = 10, path = oplsda.path)

# opermutation
opermutation <- oplsda_result@suppLs$permMN
opermutation <- as.data.frame(opermutation)
# write.csv(opermutation,paste0(oplsda.path,"/[5] OPLS-DA permutation.csv"))

colnames(opermutation)[2:3] = c("R2Y","Q2Y")
opermutation <- opermutation[,c(2,3,7)]
opermutation2 <- reshape2::melt(opermutation,3)

R2Y.value <- opermutation2$value[opermutation2$variable == "R2Y"]
Q2Y.value <- opermutation2$value[opermutation2$variable == "Q2Y"]
XR1 = 1
YR1 = max(R2Y.value)
XR2 = mean(opermutation2$sim[opermutation2$variable == "R2Y"])
YR2 = mean(R2Y.value)
slop1 <- (YR1-YR2)/(XR1-XR2)
Intecep1 <- YR1-slop1*XR1

XQ1 = 1
YQ1 = max(Q2Y.value)
XQ2 = mean(opermutation2$sim[opermutation2$variable == "Q2Y"])
YQ2 =mean(Q2Y.value)
slop2 <- (YQ1-YQ2)/(XQ1-XQ2)
Intecep2 <- YQ1-slop2*XQ1

OPLSDA_permutation <- ggscatter(opermutation2,x= "sim", y = "value",
                                fill = "variable",
                                xlab = "Permutation correlation",
                                ylab = "R2 and Q2",
                                shape = 21,
                                size = 5,
                                palette = "jco",
                                alpha = 0.8,
                                legend.title=" ",
                                font.main = c(16, "plain", "black"),            
                                font.x = c(16, "plain", "black"),                  
                                font.y = c(16, "plain", "black"),                
                                font.legend = c(16, "plain", "black")
) +
  geom_hline(yintercept = max(opermutation$R2Y), linetype = "dashed") +
  geom_hline(yintercept = max(opermutation$Q2Y), linetype = "dashed") +
  geom_abline(intercept = Intecep1, slope = slop1, color="#0077C2",
              size=0.5) +
  geom_abline(intercept = Intecep2, slope = slop2, color="#EFC000",
              size=0.5)+
  theme(axis.text = element_text(size = 16))

ggsave("[5] OPLS-DA permutation.png",
       egg::set_panel_size(OPLSDA_permutation, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = oplsda.path)
ggsave("[5] OPLS-DA permutation.pdf",
       egg::set_panel_size(OPLSDA_permutation, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5,path = oplsda.path)


# Splot
s <- as.matrix(PCA_file$marix_pareto)
T <- as.matrix(oplsda_result@scoreMN)
p1 <- NULL
for (i in 1:ncol(s)) {
  p1[i] <- as.matrix(cov(s[,i], T),ncol =1)
}

pcorr1 <-  NULL
for (i in 1:length(p1)) {
  den <- apply(T, 2, sd)*sd(s[,i])
  pcorr1[i] <- p1[i]/den
}

Splotdata <- data.frame(p1 = p1,
                        pcorr1 =pcorr1,
                        name = compound_name)

Splotdata %>% mutate(distance = p1^2+pcorr1^2) %>% 
  arrange(desc(distance)) -> Splot_distance

Splot_distance <- mutate(Splot_distance, label = name)
Splot_distance$label[11:nrow(Splot_distance)] <- NA

Splot_distance$color_label <- abs(Splot_distance$pcorr1)
Splot_distance$color_label2 <- Splot_distance$pcorr1
Splot_distance$size_label <- abs(Splot_distance$p1)

# write.csv(opermutation,paste0(oplsda.path,"/[6] OPLS-DA S-plot.csv"))

OPLSDA_Splot <- ggscatter(Splot_distance,x= "p1", y = "pcorr1",
                          fill = "color_label",
                          shape = 21,
                          legend = "top",
                          size = "size_label",
                          alpha = 0.5,
                          xlab = "p[1]",
                          ylab = "p(corr)[1]",
                          font.main = c(16, "plain", "black"),
                          font.x = c(16, "plain", "black"),
                          font.y = c(16, "plain", "black"),
                          font.legend = c(16, "plain", "black"))+
  scale_fill_gradient(low = "grey", high = "#0077C2") +
  geom_text_repel(data = Splot_distance, size = 3, aes(label = label)) +
  guides(fill = F, size = F) +
  theme(axis.text = element_text(size = 16)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("[6] OPLS-DA S-plot.png",
       egg::set_panel_size(OPLSDA_Splot, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = oplsda.path)
ggsave("[6] OPLS-DA S-plot.pdf",
       egg::set_panel_size(OPLSDA_Splot, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = oplsda.path)

OPLSDA_Splot2 <- ggscatter(Splot_distance,x= "p1", y = "pcorr1",
                          fill = "color_label2",
                          shape = 21,
                          legend = "top",
                          size = "size_label",
                          alpha = 0.5,
                          xlab = "p[1]",
                          ylab = "p(corr)[1]",
                          font.main = c(16, "plain", "black"),
                          font.x = c(16, "plain", "black"),
                          font.y = c(16, "plain", "black"),
                          font.legend = c(16, "plain", "black"))+
  scale_fill_gradient2(low = "#0077C2", mid = "white", high = "#EFC000", midpoint = 0) +
  geom_text_repel(data = Splot_distance, size = 3, aes(label = label)) +
  guides(fill = F, size = F) +
  theme(axis.text = element_text(size = 16)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("[6] OPLS-DA S-plot2.png",
       egg::set_panel_size(OPLSDA_Splot2, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = oplsda.path)
ggsave("[6] OPLS-DA S-plot2.pdf",
       egg::set_panel_size(OPLSDA_Splot2, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = oplsda.path)

wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
addWorksheet(wb, "Sheet 3", tabColour = "#008000")
addWorksheet(wb, "Sheet 4", tabColour = "#87a96b")
addWorksheet(wb, "Sheet 5", tabColour = "#ff2052")
addWorksheet(wb, "Sheet 6", tabColour = "#007fff")
addWorksheet(wb, "Sheet 7", tabColour = "#6699cc")
addWorksheet(wb, "Sheet 8", tabColour = "#de5d83")
writeData(wb, sheet = 1, oplsda_modelscore)
writeData(wb, sheet = 2, oplsda_plot)
writeData(wb, sheet = 3, oplsda_scoreMN)
writeData(wb, sheet = 4, oplsda_loadingMN_label)
writeData(wb, sheet = 5, oplsda_VIP)
writeData(wb, sheet = 6, opermutation)
writeData(wb, sheet = 7, opermutation2)
writeData(wb, sheet = 8, Splot_distance)
names(wb)[[1]] <- "oplsda modelscore"
names(wb)[[2]] <- "oplsda modelscoreplot"
names(wb)[[3]] <- "oplsda score"
names(wb)[[4]] <- "oplsda loading"
names(wb)[[5]] <- "oplsda VIP"
names(wb)[[6]] <- "oplsda opermutation"
names(wb)[[7]] <- "oplsda opermutationplot"
names(wb)[[8]] <- "oplsda Splot"

saveWorkbook(wb, paste0(oplsda.path,"/oplsda data.xlsx"), overwrite = TRUE)