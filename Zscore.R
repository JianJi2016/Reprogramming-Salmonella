dir.create(paste0(pathload,"/Zscore"))
zscore.path <- paste0(pathload,"/Zscore")

colnames(metabolites_select)
metabolites_pvalue_FC <- metabolites_select %>%
  filter(VIP > 1) %>%
  filter(ttest_p < 0.05) %>%
  filter(FC < 1/1.5 | FC > 1.5 ) %>% 
  select(name) %>%
  pull

if(length(metabolites_pvalue_FC) > 20){
  metabolites_pvalue_FC = metabolites_pvalue_FC[1:20]
}

nheigh = length(metabolites_pvalue_FC)
data_volcano_sig <- data_gapfilled %>% 
  select(c("Sample","Treatment",metabolites_pvalue_FC))

zscoreF <-function(x) {(x - mean(x))/sd(x)}
data_zscore <- apply(data_volcano_sig[,-c(1,2)],2,zscoreF)
data_zscored <- cbind(data.frame(Treatment = data_volcano_sig$Treatment),
                      data_zscore)

c <-NULL
for(i in 2:ncol(data_zscored)){
  a <-NULL
  a <- matrix(c(data_zscored[,i], rep(colnames(data_zscored)[i], nrow(data_zscored))), nrow = nrow(data_zscored), ncol = 2)
  b <- cbind(data_zscored[,1],a)
  rownames(b) = rownames(data_zscored)
  c <-rbind(c, b)
}
colnames(c) = c("Treatment", "Z_score","Metabolites")
data_zscoreReady = c

for (i in 1:length(levels(data_volcano_sig$Treatment))){
  print(paste0(i))
  data_zscoreReady[,1][data_zscoreReady[,1] == paste0(i)] = levels(data_volcano_sig$Treatment)[i]
}
data_zscoreReady = data.frame(data_zscoreReady, stringsAsFactors = F)

data_zscoreReady$Z_score = as.numeric(data_zscoreReady$Z_score)

zscore <- data_zscoreReady %>%
  ggscatter(.,
            x = "Metabolites", 
            y = "Z_score",
            size = 5,
            xlab = "",
            ylab = "Z-score",
            ylim = c(-3, 3),
            # combine = T,
            shape = 16,
            color = "Treatment", 
            palette = "jco", 
            # rotate = T,
            font.main = c(16, "plain", "black"), 
            font.x = c(16, "plain", "black"),                  
            font.y = c(16, "plain", "black"),                
            legend = "top", 
            font.legend = c(16, "plain", "black"), 
            repel = T) +
  coord_flip() +
  geom_vline(xintercept = c(a <- c(1:nrow(data_zscoreReady))), 
             linetype = "dashed") +
  theme(axis.text = element_text(size = 16)) +
  guides(fill = guide_legend(title = ''))

ggsave("[1] Z-score plsda VIP10.png", 
       egg::set_panel_size(zscore, 
                           width=unit(3.5, "in"), 
                           height=unit(6, "in")), 
       height = 8, 
       width = 10,path = zscore.path)

ggsave("[1] Z-score plsda VIP10.pdf", 
       egg::set_panel_size(zscore, 
                           width=unit(3.5, "in"), 
                           height=unit(6, "in")), 
       height = 8, 
       width = 10,path = zscore.path)


wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
addWorksheet(wb, "Sheet 3", tabColour = "#008000")
writeData(wb, sheet = 1, data_volcano_sig)
writeData(wb, sheet = 2, data_zscored)
writeData(wb, sheet = 3, data_zscoreReady)
names(wb)[[1]] <- "data_volcano_sig"
names(wb)[[2]] <- "data_zscored"
names(wb)[[3]] <- "data_zscoreReady"
saveWorkbook(wb, paste0(zscore.path,"/data Zscore.xlsx"), overwrite = TRUE)
