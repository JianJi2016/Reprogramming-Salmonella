dir.create(paste0(pathload,"/Heatmap"))
Heatmap.path <- paste0(pathload,"/Heatmap")

annotation_col <- data_volcano_sig %>% 
  select("Sample","Treatment")
rownames(annotation_col) = data_volcano_sig[,1]
annotation_col = annotation_col %>% select(Treatment)

test <- data_volcano_sig %>% 
  select(-c("Sample","Treatment"))
rownames(test) = data_volcano_sig[,1]

ann_colors = list(Treatment = c(L = '#0077C2',
                          R = '#EFC000')) 
names(ann_colors$Treatment) = c(data_volcano_sig$Treatment[1],
                            data_volcano_sig$Treatment[nrow(data_volcano_sig)])

test2 <- test %>%
  log10()
colnames(test2) <- colnames(test2) %>%
  substr(.,1,30)

test_b <- pheatmap(t(test2),
                   scale = "row",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   annotation_col = annotation_col,
                   # annotation_row = annotation_row,
                   annotation_colors = ann_colors,
                   # labels_row = labels_row,
                   border=FALSE,
                   border_color = "black",
                   treeheight_row = 30, 
                   treeheight_col = 30,
                   cellwidth = 12,
                   cellheight = 12,
                   fontsize = 8,
                   clustering_method = "ward.D",
                   # "ward.D","single", "complete", "ward.D2"
                   # cluster_rows = FALSE, 
                   cluster_cols = F
                   # cutree_col = 2
                   # display_numbers = TRUE,
                   # display_numbers = matrix(ifelse(test > 0.05, "*", ""), nrow(test))
)
pheatmap = ggplotify::as.ggplot(test_b)

nrow(test)
ncol(test)
ggsave("pheatmap2.pdf", 
       egg::set_panel_size(pheatmap, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       width = 20, height = 5,
       path = Heatmap.path)

ggsave("pheatmap2.png", 
       egg::set_panel_size(pheatmap, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       width = 20, height = 5,
       path = Heatmap.path)


# -------------


metabolites_select %>%
  arrange(ttest_p) %>% top_n(., n = 20) %>% select(name) %>% pull -> forheatmap

annotation_col <- data_volcano_sig %>% 
  select("Sample","Treatment")
rownames(annotation_col) = data_volcano_sig[,1]
annotation_col = annotation_col %>% select(Treatment)

test <- data_gapfilled %>% 
  select(forheatmap)
rownames(test) = data_volcano_sig[,1]

ann_colors = list(Treatment = c(L = '#0077C2',
                                R = '#EFC000')) 
names(ann_colors$Treatment) = c(data_volcano_sig$Treatment[1],
                                data_volcano_sig$Treatment[nrow(data_volcano_sig)])
test2 <- test %>%
  log10()
colnames(test2) <- colnames(test2) %>%
  substr(.,1,30)
test_b <- pheatmap(t(test2),
                   scale = "row",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   annotation_col = annotation_col,
                   # annotation_row = annotation_row,
                   annotation_colors = ann_colors,
                   # labels_row = labels_row,
                   border=FALSE,
                   border_color = "black",
                   treeheight_row = 30,
                   treeheight_col = 30,
                   cellwidth = 12,
                   cellheight = 12,
                   fontsize = 8,
                   clustering_method = "ward.D",
                   cluster_cols = FALSE
)
pheatmap = ggplotify::as.ggplot(test_b)
ggsave("pheatmap.pdf", 
       egg::set_panel_size(pheatmap, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       width = 20, height = 5,
       path = Heatmap.path)

ggsave("pheatmap.png", 
       egg::set_panel_size(pheatmap, 
                           width=unit(3.5, "in"), 
                           height=unit(3.5, "in")), 
       width = 20, height = 5,
       path = Heatmap.path)

# ------
wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
writeData(wb, sheet = 1, test)
writeData(wb, sheet = 2, annotation_col)
names(wb)[[1]] <- "test"
names(wb)[[2]] <- "annotation_col"
saveWorkbook(wb, paste0(Heatmap.path,"/data Heatmap.xlsx"), overwrite = TRUE)
