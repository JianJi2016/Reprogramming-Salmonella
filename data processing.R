cat("import raw csv data \n")
file <- list.files(pathload)
sample_raw <- read.csv(paste0(pathload,"/",file),check.names = F)

sample_import <- sample_raw


cat("transfer all peak intensity to numeric format \n")
sample_import[,3:ncol(sample_import)] <- sample_import[,3:ncol(sample_import)] %>%
  lapply(., as.numeric) %>%
  as.data.frame()

cat("change NA to 0 \n")
sample_import[is.na(sample_import)] <- 0

# change all the value to positive
sample_import[,3:ncol(sample_import)] <-  abs(sample_import[,3:ncol(sample_import)])

# removal the metabolites with over 30% zero value
zerocheck <- NULL
zerocheck[1:2] = c(0,0)
for (i in 3:ncol(sample_import)) {
  sample_import[sample_import[,i] == 0, i] %>% length() ->  zerovalue
  sample_import[, i] %>% length() -> all
  zerocheck[i] = zerovalue/all
}
sample_import = sample_import[, zerocheck < 0.3]

# check the metabolites with same value
samecheck <- NULL
samecheck[1:ncol(sample_import)] = rep(T,ncol(sample_import))
for (i in 3:ncol(sample_import)) {
  if(max(sample_import[,i]) == min(sample_import[,i])){
    samecheck[i] = F
  }
}
sample_import = sample_import[, samecheck]

# rename the first two column and format the treatment group
colnames(sample_import)[1:2] <- c("Sample","Treatment")



sample_import$Treatment <- factor(sample_import$Treatment,
                                  levels = c(
                                    sample_import$Treatment[1],
                                    sample_import$Treatment[nrow(sample_import)]
                                            )
                                  )


# data exported after cleaning
# write.csv(sample_import, "[1] data_clean.csv", row.names = F)

sampleinfo <- sample_import %>% 
  select(.,Sample, Treatment) %>%
  textshape::column_to_rownames(.)

colnames(sample_import)[-c(1:2)] %>% 
  as.data.frame() -> medi
colnames(medi) = "Compounds"
metainfo <- medi

compound_name <- colnames(sample_import)[-c(1:2)]

matrix <- sample_import %>% 
  select(.,-Sample, -Treatment)

# write.csv(sampleinfo, "[2] sampleinfo.csv", row.names = F)
# write.csv(metainfo, "[3] metainfo.csv", row.names = F)
# write.csv(matrix, "[4] matrix.csv", row.names = F)

#  gap fill, use the min value of the metabolites in all samples
matrix_fill <- matrix
for (i in 1:ncol(matrix_fill)) {
  if(grep("TRUE", matrix_fill[,i] == 0) %>% length() > 0){
    matrix_fill[,i][grep("TRUE", matrix_fill[,i] == 0)] <-
      matrix_fill[,i][grep("TRUE", matrix_fill[,i] != 0)] %>% min()
  }
}

data_gapfilled <- cbind(sample_import[,c(1:2)],matrix_fill)
# write.csv(data_gapfilled, "[5] data_gapfilled.csv", row.names = F)

# log transformation
matrix_log <- matrix_fill %>% log10
data_log <- matrix_log %>%
  cbind(sample_import[,c(1:2)],.)
# write.csv(data_log, "[6] data_log.csv", row.names = F)

# scaling data
pareto_scale <- function(x){(x-mean(x, na.rm = T))/sqrt(sd(x, na.rm = T))}
# auto_scale <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

marix_pareto  <- apply(matrix_log, 2, pareto_scale)
data_pareto <- cbind(sample_import[,c(1:2)], marix_pareto)
# write.csv(data_pareto, "[7] data_pareto.csv", row.names = F)

hs <- createStyle(
  textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,
  fontName = "Arial Narrow", fgFill = "#4F80BD"
)

wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
addWorksheet(wb, "Sheet 2", tabColour = "#0070ff")
addWorksheet(wb, "Sheet 3", tabColour = "#008000")
addWorksheet(wb, "Sheet 4", tabColour = "#87a96b")
addWorksheet(wb, "Sheet 5", tabColour = "#ff2052")
addWorksheet(wb, "Sheet 6", tabColour = "#007fff")
addWorksheet(wb, "Sheet 7", tabColour = "#6699cc")
addWorksheet(wb, "Sheet 8", tabColour = "#de5d83")
writeData(wb, sheet = 1, sample_raw)
writeData(wb, sheet = 2, sample_import)
writeData(wb, sheet = 3, sampleinfo)
writeData(wb, sheet = 4, metainfo)
writeData(wb, sheet = 5, matrix)
writeData(wb, sheet = 6, data_gapfilled)
writeData(wb, sheet = 7, data_log)
writeData(wb, sheet = 8, data_pareto)
names(wb)[[1]] <- "sample_raw"
names(wb)[[2]] <- "sample_import"
names(wb)[[3]] <- "sampleinfo"
names(wb)[[4]] <- "metainfo"
names(wb)[[5]] <- "matrix"
names(wb)[[6]] <- "data_gapfilled"
names(wb)[[7]] <- "data_log"
names(wb)[[8]] <- "data_pareto"

dir.create(paste0(pathload,"/data processed"))
process.path <- paste0(pathload,"/data processed")
saveWorkbook(wb, file = paste0(process.path,"/data processed.xlsx"), overwrite = TRUE)


