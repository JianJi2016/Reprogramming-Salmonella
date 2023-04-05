dir.create(paste0(pathload,"/Volcano Results"))
volcano.path <- paste0(pathload,"/Volcano Results")

data_gapfilled$Treatment <- factor(data_gapfilled$Treatment,
                                   levels = data_gapfilled %>% 
                                     select(Treatment) %>% 
                                     unique() %>% 
                                     pull(1)
)

factorNumber <- length(levels(data_gapfilled$Treatment))
data_original <- data_gapfilled

colnames(data_original)[-c(1,2)] <- compound_name
factors <- c(data_original$Treatment[1],data_original$Treatment[nrow(data_original)])
treatments <- data_original$Treatment %>% unique() %>% length()
TreatmentA <- which(data_original$Treatment == factors[1]) ## non-cancer tissues
TreatmentB <- which(data_original$Treatment == factors[2]) ## tumor tissues
runttest <- function(data1, data2) {
  results <- t.test(data1,data2)
  results2 <-  wilcox.test(data1,data2)
  foldchange <- mean(data2)/mean(data1)
  unlist(c(results$p.value,results2$p.value,foldchange))
}

ttest_value <- NULL
utest_value <- NULL
FC_value <- NULL
for(i in 3:ncol(data_original)){
  ttest_value[i-2] = runttest(data_original[TreatmentA,i],data_original[TreatmentB,i])[1]
  utest_value[i-2] = runttest(data_original[TreatmentA,i],data_original[TreatmentB,i])[2]
  FC_value[i-2] = runttest(data_original[TreatmentA,i],data_original[TreatmentB,i])[3]
}

data_original_tt_FC_result <- data.frame(Sample =colnames(data_original)[-c(1,2)], 
                                         pvalue = ttest_value,
                                         pvalue2 = utest_value,
                                         FC = FC_value, 
                                         p.adjust = p.adjust(ttest_value, method = "fdr"),
                                         p.adjust2 = p.adjust(utest_value, method = "bonferroni")
)
vacanno <- data_original_tt_FC_result
colnames(vacanno) = c("compound","ttest_p","utest_p","FC","fdr_adjust_ttest_p",
                      "bonferroni_adjust_utest_p")
vacanno <- vacanno %>%
  mutate(log2FoldChange = log2(FC)) %>%
  mutate(pvalue = ttest_p)

de <- vacanno
# ttest for plot
value = 1.5
de$diffexpressed <- "NO"
de$diffexpressed[de$FC > value & de$pvalue < 0.05] <- "UP"
de$diffexpressed[de$FC < 1/value & de$pvalue < 0.05] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$compound[de$diffexpressed != "NO"]

# de$log2FoldChange <- de$FoldChange %>% log(.,2)
# de$Logpvalue2 <- de$pvalue %>% log10()
# de$Logpvalue <- de$Logpvalue2 * -1

volcanoplot <- ggplot(data = de,
                      size = 5,
                      aes(x = log2FoldChange, 
                      y = -log10(pvalue), 
                      fill = diffexpressed,
                      label= delabel)) + 
  geom_text_repel() +
  geom_point(shape = 21, size = 3) + ggpubr::theme_pubr() +
  scale_fill_manual(values =c("UP" = "#0077C2", 
                               "NO" = "grey", 
                               "DOWN" = "#EFC000")) + 
  geom_vline(xintercept = c(log2(1/value), log2(value)), col = "grey", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "longdash") +
  labs(fill = "")

ggsave("[1] volcano_ttest.png",
       egg::set_panel_size(volcanoplot, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = volcano.path)

ggsave("[1] volcano_ttest.pdf",
       egg::set_panel_size(volcanoplot, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, width = 5, path = volcano.path)

# utest for plot
de$diffexpressed2 <- "NO"
de$diffexpressed2[de$FC > value & de$utest_p < 0.05] <- "UP"
de$diffexpressed2[de$FC < 1/value & de$utest_p < 0.05] <- "DOWN"
de$delabel2 <- NA
de$delabel2[de$diffexpressed2 != "NO"] <- de$compound[de$diffexpressed2 != "NO"]

volcanoplot2 <- ggplot(data = de, aes(x = log2FoldChange, 
                      y = -log10(utest_p), 
                      fill = diffexpressed2,
                      label= delabel2)) + 
  geom_text_repel() +
  geom_point(shape = 21, size = 3) + ggpubr::theme_pubr() +
  scale_fill_manual(values =c("UP" = "#0077C2", 
                              "NO" = "grey", 
                               "DOWN" = "#EFC000")) + 
  geom_vline(xintercept = c(log2(1/value), log2(value)), col = "grey", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "longdash") +
  labs(fill = "")

ggsave("[2] volcano_utest.png",
       egg::set_panel_size(volcanoplot2, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height =5, width = 5, path = volcano.path)

ggsave("[2] volcano_utest.pdf",
       egg::set_panel_size(volcanoplot2, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height =5, width = 5, path = volcano.path)




#  barplot
vacanno5 <- de %>%
  filter(diffexpressed != "NO") %>%
  mutate(ABSlog2FoldChange = log2FoldChange %>% abs()) %>%
  arrange(desc(ABSlog2FoldChange)) %>%
  top_n(n = 20)

volbar <- ggbarplot(vacanno5 %>%
          filter(diffexpressed != "NO"), 
          x = "compound", y = "log2FoldChange",
          color = "black",            # Set bar border colors to white
          sort.val = "asc",           # Sort the value in ascending order
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Log2(FC value)",
          xlab = " ",
          rotate = TRUE,
          fill = "log2FoldChange", 
          legend.title = " ",
          legend = " ",
          font.label = list(color = "white", size = 9, 
                            vjust = 0)) +
  scale_fill_gradient2(low = "#0077C2", 
                       mid = "white", 
                       high = "#EFC000", 
                       midpoint = 0) +
  geom_hline(yintercept = c(a <- c(log2(value),log2(1/value))), 
             linetype = "dashed")  

ggsave("[3] Volcano Barplot.png", 
       egg::set_panel_size(volbar, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, 
       width = 10, path = volcano.path)
ggsave("[3] Volcano Barplot.pdf", 
       egg::set_panel_size(volbar, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, 
       width = 10, path = volcano.path)


nUP <- de %>%
  filter(diffexpressed == "UP") %>%
  nrow()

nDOWN <- de %>%
  filter(diffexpressed == "DOWN") %>%
  nrow()

nALL <- de %>% nrow()

sigdata <- data.frame(up = nUP,
           down = nDOWN,
           all = nALL)

sigdata <- sigdata %>%
  mutate(none = all-up-down)
sigdataF <- t(sigdata)
colnames(sigdataF) = "number"
sigdataF <- sigdataF %>%
  as.data.frame() %>%
  rownames_to_column() 


Bartren <- ggbarplot(sigdataF %>%
            filter(rowname != "all"),
          x = "rowname",
          y = "number",
          fill = "rowname",
          palette = "jco",
          label = TRUE, label.pos = "out"
          )



ggsave("[4] Volcano Bartren.png", 
       egg::set_panel_size(Bartren, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, 
       width = 5, path = volcano.path)
ggsave("[4] Volcano Bartren.pdf", 
       egg::set_panel_size(Bartren, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       height = 5, 
       width = 5, path = volcano.path)



wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
addWorksheet(wb, "Sheet 3", tabColour = "#008000")
addWorksheet(wb, "Sheet 4", tabColour = "#008000")

writeData(wb, sheet = 1, data_original)
writeData(wb, sheet = 2, de)
writeData(wb, sheet = 3, vacanno5)
writeData(wb, sheet = 4, sigdata)

names(wb)[[1]] <- "raw data"
names(wb)[[2]] <- "scatter plot"
names(wb)[[3]] <- "barplot"
names(wb)[[4]] <- "sigdata"
saveWorkbook(wb, paste0(volcano.path,"/volcano data.xlsx"), overwrite = TRUE)

