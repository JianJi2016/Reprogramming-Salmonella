cat("---------------------|Loading packages|--------------------------")
cat("\n") 

pacman::p_load(textshape,ropls,ggpubr,ggplot2,ggrepel,tibble,stats,
               ggsignif,pROC,plotROC,
               circlize,ggcorrplot,stringr,tidyr,factoextra,
               readr,openxlsx,dplyr,pheatmap,ggplotify,egg,corrplot)
rm(list=ls())
options(warn = -1)

pathload <- getwd()
setwd("F:/MetaPro ver 3.0")

cat("---------------------|read data and data transformation|---")
cat("\n") 
source("data processing.R")

cat("---------------------|PCA starting|------------------------")
cat("\n") 
source("PCA.R")


cat("---------------------|OPLS-DA starting|--------------------")
cat("\n") 
source("OPLS-DA.R")

cat("---------------------|volcano starting|--------------------")
cat("\n") 
source("Volcano.R")

cat("---------------------|Biomarkers|--------------------------")
cat("\n") 
source("Biomarkers.R")

cat("---------------------|Boxplot starting|--------------------")
cat("\n")
source("Boxplot.R")

cat("---------------------|Zscore starting|---------------------")
cat("\n") 
source("Zscore.R")

cat("---------------------|ROC starting|--------------------------")
cat("\n")
source("ROC.R")

cat("---------------------|Heatmap starting|--------------------------")
cat("\n") 
source("Heatmap.R")


cat("---------------------|Corrplot starting|--------------------------")
cat("\n") 
source("Corrplot.R")

cat("---------------------|clustering ploting|--------------------------")
cat("\n") 
source("Cluster.R")

cat("---------------------|Samplecor|--------------------------")
cat("\n") 
source("Samplecor.R")

cat("---------------------|Pathway enrichment starting|--------------------------")
cat("\n")
source("Pathway.R")

