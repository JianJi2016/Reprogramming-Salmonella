dir.create(paste0(pathload,"/Boxplot"))
boxplot.path <- paste0(pathload,"/Boxplot")

metabolites_selectinfo <- vacanno %>% 
  filter(ttest_p < 0.05 | utest_p < 0.05) 

compound_select <- vacanno %>% 
  filter(ttest_p < 0.05 | utest_p < 0.05) %>%
  select(compound)  %>%
  pull()

data_select <- data_gapfilled %>%
  select(c("Sample","Treatment",compound_select))

colnames(data_select) = colnames(data_select) %>% make.names()
compound_select = make.names(compound_select)

for (i in 1:length(compound_select)) {
  cat("making boxplot, No: ",i,"\n")

  boxplotN <- ggpar(ggboxplot(data_select, 
                  x = "Treatment", 
                  y = compound_select[i+2], 
                  color = "Treatment",
                  # fill = "Treatment",
                  add = c("point"),
                  xlab = "Treatment",
                  width = 0.6,
                  ylab = "Relative Intensity",
                  add.params = list(size = 3, 
                                    jitter = 0.2)),
        font.main = c(16, "plain", "black"),            
        font.x = c(16, "plain", "black"),                  
        font.y = c(16, "plain", "black"),                
        legend = " ",                                  
        font.legend = c(16, "plain", "black"),               
        rotate = F,                                         
        # ticks = T,                                          
        # tickslab = T, 
        palette = "jco",
        xtickslab.rt = 0) +
    theme(axis.text = element_text(size = 16)) +
    stat_signif(comparisons = list(data_select$Treatment %>% levels()),
                test = "wilcox.test")
  
  ggsave(paste0(make.names(compound_select[i]), ".png"),
         egg::set_panel_size(boxplotN, 
                             width=unit(3.5, "in"), 
                             height=unit(3.5, "in")), 
         height = 5, width = 5,path = boxplot.path)
  
  
  ggsave(paste0(make.names(compound_select[i]), ".pdf"),
         egg::set_panel_size(boxplotN, 
                             width=unit(3.5, "in"), 
                             height=unit(3.5, "in")), 
         height = 5, width = 5,path = boxplot.path)
  
}

hs <- createStyle(
  textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,
  fontName = "Arial Narrow", fgFill = "#4F80BD"
)

wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
addWorksheet(wb, "Sheet 3", tabColour = "#008000")
writeData(wb, sheet = 1, compound_select)
writeData(wb, sheet = 2, data_select)
writeData(wb, sheet = 3, metabolites_selectinfo)
names(wb)[[1]] <- "compound_select"
names(wb)[[2]] <- "data_select"
names(wb)[[3]] <- "metabolites_selectinfo"
saveWorkbook(wb, paste0(boxplot.path,"/data boxplot.xlsx"), overwrite = TRUE)