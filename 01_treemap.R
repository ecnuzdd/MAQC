#### set work directory ####
BASE_DIR = normalizePath(
  "D:/MAQC/20200508/Analysis/code_library", mustWork = F
)
setwd(BASE_DIR)

# library
library(xlsx)
library(readxl)
library(stringr)
library(treemap)

# input
tm_df.path = normalizePath(
  file.path(BASE_DIR, 'temp_figures', '01', 'analysis_CC_by_tiansha.xlsx'), mustWork = F
)
tm_df = read_excel(tm_df.path, sheet = 'treemap')
tm_df$logFDR = -log10(tm_df$FDR)
tm_df$CC1 = paste(
  tm_df$CC, ' (', tm_df$GPs, ')', sep = ''
)

tm_df.plot <- treemap(
  tm_df,
  index=c("CC1"),
  vSize="GPs",
  type="index",
  palette = "Pastel1",
  bg.labels=c("white"),
  align.labels=list(
    c("center", "center")
  ),
  title = 'Cellular component of Quartet\n(The union set: 12,068 GPs)'
)            


library(treemap)
df = data.frame(
  Comment = c(
    'ABC transportor', 'Chromatin', 'DNA replication', 'Enzyme',
    'RNA', 'Subcellular localization', 'Trafficking', 'Others'
  ),
  GPs = c(3, 6, 9, 13, 17, 19, 3, 3)
)
df$Comment = paste(
  df$Comment, ' (', df$GPs, ')', sep = ''
)
tm_df.plot <- treemap(
  df,
  index=c("Comment"),
  vSize="GPs",
  type="index",
  palette = "Pastel1",
  bg.labels=c("white"),
  align.labels=list(
    c("center", "center")
  ),
  title = 'CSG-overexpressing proteins (Total = 71)'
)    

