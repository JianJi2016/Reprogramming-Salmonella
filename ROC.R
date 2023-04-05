dir.create(paste0(pathload,"/ROC"))
roc.path <- paste0(pathload,"/ROC")

ROC_group1 <- data_volcano_sig %>% 
  filter(Treatment == data_volcano_sig$Treatment %>%
           unique() %>%
           .[1] %>% 
           as.character()
  ) 

ROC_group2 <- data_volcano_sig %>%
  filter(Treatment == data_volcano_sig$Treatment %>%
           unique() %>%
           .[2] %>% 
           as.character()
  ) 

ROC_group1$Treatment = "1"
ROC_group2$Treatment = "0"
ROC_ready <- rbind(ROC_group1,ROC_group2)

MetaROC <- NULL
for (i in 3:ncol(ROC_ready)) {
  ROCvalue <- roc(ROC_ready[,2],
                  ROC_ready[,i],
                  plot = F,
                  levels=c("0", "1"),
                  direction = ">")
  
  ROC_value <- ROCvalue$auc[1]
  MetaROC[i-2] = round(ROCvalue$auc[1],2)
  cat("Show the metabolites AUC over 0.7")
  cat("\n")
  if(ROC_value > 0.1) {
    cat(paste0(colnames(ROC_ready)[i]," AUC: ",round(ROCvalue$auc[1],2)))
    cat("\n")
    ROCplot <- pROC::ggroc(ROCvalue,
                           legacy.axes = TRUE,
                           color = "#3B4992",
                           size = 0.5) +
      labs(x = "1-Specificity", y = "Sensitivity") +
      ggpubr::theme_pubr() +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                   color = "darkgrey", linetype = "dashed") +
      annotate(geom = "text",
               x = 0.7, y = 0.1,
               size = 5,
               color = "darkgrey",
               label = paste0("AUC = ",round(ROC_value,2))) +
      theme(axis.title.y = element_text(size = 16, angle = 90)) +
      theme(axis.title.x = element_text(size = 16, angle = 0)) +
      theme(axis.text = element_text(size = 16))
    print(ROCplot)
    ggsave(paste0(make.names(colnames(ROC_ready)[i]),".png"),
           egg::set_panel_size(ROCplot, 
                               width=unit(3.5, "in"), 
                               height=unit(3.5, "in")), 
           height = 5, width = 5, path =roc.path )
    pdf(paste0(roc.path,"/",make.names(colnames(ROC_ready)[i]),".pdf"),
        egg::set_panel_size(ROCplot, 
                            width=unit(3.5, "in"), 
                            height=unit(3.5, "in")), 
        height = 5, width = 5)
    print(ROCplot)
    dev.off()
  }
}

MetaROC <- as.data.frame(MetaROC)
rownames(MetaROC) = colnames(ROC_ready)[-(1:2)]
colnames(MetaROC)[1] = "AUC"

MetaROC$order <- c(1:nrow(MetaROC))
MetaROC$name <- rownames(MetaROC)
MetaROC$label <- ""
MetaROC$fill <- " "
for(i in 1:nrow(MetaROC)){
  if(MetaROC$AUC[i] > 0.9){
    MetaROC$label[i] = MetaROC$name[i]
    MetaROC$fill[i] = "red"
  }else{
    MetaROC$label[i] = NA
    MetaROC$fill[i] = "grey"
  }
}
MetaROC$fill <- as.factor(MetaROC$fill)

ROC_plot <- ggscatter(MetaROC,"order","AUC",
                      shape = 21,fill = "fill",
                      color = "black",size = "AUC",
                      palette = "aaas",
                      alpha = "fill",
                      xlab = "Metabolites",
                      ylab = "AUC score",
                      font.main = c(16, "plain", "black"), 
                      font.x = c(16, "plain", "black"),                  
                      font.y = c(16, "plain", "black"),                
                      # legend = "top", 
                      # font.legend = c(16, "plain", "black"),
                      legend = "") +
  geom_text_repel(
    data = MetaROC,
    aes(x=order,y=AUC,label = label),
    size = 3,
    segment.color = "black", show.legend = FALSE ) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  theme(axis.text = element_text(size = 16))
ggsave("[1] AUC plot.png",height = 5, width = 5, path =roc.path)
ggsave("[1] AUC plot.pdf",height = 5, width = 5, path =roc.path)

MetaROC %>% arrange(AUC) -> MetaROC

wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
writeData(wb, sheet = 1, MetaROC)
writeData(wb, sheet = 2, ROC_ready)
names(wb)[[1]] <- "MetaROC"
names(wb)[[2]] <- "ROC_ready"
saveWorkbook(wb, paste0(roc.path,"/data ROC.xlsx"), overwrite = TRUE)
