dir.create(paste0(pathload,"/Pathway enrichment"))
Pathway.path <- paste0(pathload,"/Pathway enrichment")

file.list <- dir(Pathway.path)
length(file.list)
if(length(file.list) == 1){
  pathway.data <- read_csv(paste0(Pathway.path,"/",file.list))
  pacman::p_load(ggpubr, ggplot2, dplyr)
  
  colnames(pathway.data)[1] <- "Name"
  colnames(pathway.data)[6] <- "logP"
  
  Position<- order(pathway.data$Impact,decreasing=TRUE)[1:3]
  pathway.data$Label <- " "
  pathway.data$Label[Position] <- as.character(pathway.data$Name[Position])
  
  PathwayFig <- ggscatter(pathway.data, x = "Impact", y = "logP",
                          size = "Impact", fill = "logP",
                          shape = 21,color = "black",
                          font.label = 8, repel = TRUE,
                          xlab = "Pathway Impact",
                          ylab = "-ln P-value",
                          legend = "right",
                          legend.title = "-ln P-value",
                          label = "Label",
                          font.x = c(14, "plain", "black"),
                          font.y = c(14,  "plain", "black"),
                          font.legend = c(14, "plain", "black"))+
    scale_fill_gradient(low = "yellow", high = "red") +
    theme(legend.position="right") +
    labs(size = "Pathway Impact") +
    guides(size = guide_legend(order = 1))
  
  ggsave("[1] PathwayFig.png",height = 4, width = 5, path = Pathway.path)
  ggsave("[1] PathwayFig.pdf",height = 4, width = 5, path = Pathway.path)
  
  PathwayBar <- pathway.data %>%
    head(15) %>%
    arrange(Hits) %>%
    ggbarplot(.,x = "Name", y = "Hits",
              xlab = " ",
              fill = "Hits",
              ylab = "Number of metabolites") +
    coord_flip()
  
  ggsave("[2] PathwayBar.png",height = 6, width = 6, path = Pathway.path)
  ggsave("[2] PathwayBar.pdf",height = 6, width = 6, path = Pathway.path)
  
  
  PathwayDot <- pathway.data %>%
    head(12) %>%
    arrange(Hits) %>%
    ggscatter(.,x = "Name", y = "Impact",
              size = "Hits",
              color = "logP",
              xlab = " ",
              ylab = "Impact factor") +
    coord_flip() +
    geom_vline(xintercept = c(a <- c(1:nrow(pathway.data))), linetype = "dashed") +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_binned(range = c(2, 10))
  
  ggsave("[3] PathwayDot.png",height = 6, width = 6, path = Pathway.path)
  ggsave("[3] PathwayDot.pdf",height = 6, width = 6, path = Pathway.path)
}
