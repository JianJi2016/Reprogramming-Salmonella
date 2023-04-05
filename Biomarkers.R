dir.create(paste0(pathload,"/Biomarkers"))
biomarkers.path <- paste0(pathload,"/Biomarkers")

cbind(oplsda_VIP, vacanno ) -> metabolites_select
metabolites_select <- metabolites_select[,-3]

wb <- createWorkbook("Fred")
addWorksheet(wb, "Sheet 1")
writeData(wb, sheet = 1, metabolites_select)
names(wb)[[1]] <- "Metabolites"
saveWorkbook(wb, paste0(biomarkers.path,"/Metabolites.xlsx"), overwrite = TRUE)