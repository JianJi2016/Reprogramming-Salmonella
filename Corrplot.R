dir.create(paste0(pathload,"/Corrplot"))
corrplot.path <- paste0(pathload,"/Corrplot")

cname <- colnames(test) 
for (i in 1:length(cname)) {
  cname[i] = substr(cname[i],1,20)
}
colnames(test) = cname

corr <- round(cor(test), 1)
ggcorr_plot <- ggcorrplot(corr, hc.order = TRUE, 
                          type = "lower",
                          method = "square",
                          lab = TRUE,
                          colors = c("navy", "white", "firebrick3")) +
  scale_y_discrete(position = "right") +  theme_classic() +
  theme(axis.title.y = element_text(size = 12, angle = 0)) +
  theme(axis.title.x = element_text(size = 12, angle = 0)) +
  theme(axis.text.x = element_text(colour = "black", size = 12, 
                                   angle = -90, vjust = 0.5, hjust=0)) +
  theme(axis.text.y = element_text(colour = "black", size = 12, angle = 0))+
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 12))+
  labs(
    x = "",
    y = "")

Ggcorr_plot <- ggplotify::as.ggplot(ggcorr_plot)

ggsave("Corrplot_top20.pdf", 
       egg::set_panel_size(Ggcorr_plot, 
                           width=unit(nheigh/2 + 5, "in"), 
                           height=unit(nheigh/2 + 5, "in")), 
       width = 15, height = 15,
       path = corrplot.path)

ggsave("Corrplot_top20.png", 
       egg::set_panel_size(Ggcorr_plot, 
                           width=unit(nheigh/2 + 5, "in"), 
                           height=unit(nheigh/2 + 5, "in")), 
       width = 15, height = 15,
       path = corrplot.path)


wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")

names <- rownames(corr)
rownames(corr) <- NULL
corr2 <- cbind(names,corr)

writeData(wb, sheet = 1, corr2)
names(wb)[[1]] <- "Correlation"
saveWorkbook(wb, paste0(corrplot.path,"/Metabolites Correlation.xlsx"), overwrite = TRUE)
