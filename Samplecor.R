dir.create(paste0(pathload,"/Sample cor"))
Samplecor.path <- paste0(pathload,"/Sample cor")




sampleimp <- sample_import[,-2]
rownames(sampleimp) = sampleimp[,1]
sampleimp2 = sampleimp[,-1]
corr <- round(cor(t(sampleimp2)), 2)


png( 
  filename = paste0(Samplecor.path,"/samplecor.png"), # 文件名称
  width = 60,           # 宽
  height = 60,          # 高
  units = "in",          # 单位
  bg = "white",          # 背景颜色
  res = 120)              # 分辨率
corrplot(
  # 相关系数矩阵
  corr = corr, 
  order = 'AOE',
  method = 'number',
  type = 'lower', 
  tl.pos = 'tp',
  tl.col = "black")

corrplot(
  corr = corr, 
  add = TRUE, 
  type = 'upper', 
  method = 'ellipse',
  order = 'AOE',
  diag = FALSE,  
  tl.pos = 'n', 
  # cl.pos = 'n',
  tl.col = "black")

dev.off()


pdf(paste0(Samplecor.path,"/samplecor.pdf"), height = 30, width = 30)
corrplot(
  # 相关系数矩阵
  corr = corr, 
  order = 'AOE',
  method = 'number',
  type = 'lower', 
  tl.pos = 'tp',
  tl.col = "black")

corrplot(
  corr = corr, 
  add = TRUE, 
  type = 'upper', 
  method = 'ellipse',
  order = 'AOE',
  diag = FALSE,  
  tl.pos = 'n', 
  # cl.pos = 'n',
  tl.col = "black")

dev.off()

