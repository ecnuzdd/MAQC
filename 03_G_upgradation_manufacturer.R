# fot_data
# meta_data

meta_data.TMO.E480 = meta_data[which(
  meta_data$Labcode == 'TMO' 
  & meta_data$SInstrument == 'Exploris480'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]

meta_data.BRK.timsTOF = meta_data[which(
  meta_data$Labcode == 'BRK' 
  & meta_data$SInstrument == 'timsTOF'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]

meta_data.BRK.QTOF6600 = meta_data[which(
  meta_data$Labcode == 'ABS' 
  & meta_data$SInstrument == 'QTOF6600'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]




# D5
px_list = c('D5', 'D6', 'F7', 'M8')
px = px_list[1]

e480.px.index = which(meta_data.TMO.E480$UPerson == px)
e480.px.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.E480$Exp_code[e480.px.index])

timstof.px.index = which(meta_data.BRK.timsTOF$UPerson == px)
timstof.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BRK.timsTOF$Exp_code[timstof.px.index])

qtof6600.px.index = which(meta_data.BRK.QTOF6600$UPerson == px)
qtof6600.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BRK.QTOF6600$Exp_code[qtof6600.px.index])



#1 identification
e480.px.col_gps = apply(e480.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(e480.px.fot[,1])))

timstof.px.col_gps = apply(timstof.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(timstof.px.fot[,1])))

qtof6600.px.col_gps = apply(qtof6600.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(qtof6600.px.fot[,1])))


e480.px.id = as.vector(unlist(lapply(e480.px.col_gps, length)))
timstof.px.id = as.vector(unlist(lapply(timstof.px.col_gps, length)))
qtof6600.px.id = as.vector(unlist(lapply(qtof6600.px.col_gps, length)))

barplot_df = data.frame(
  ID =  c(
    e480.px.id, nrow(e480.px.fot), 
    timstof.px.id, nrow(timstof.px.fot),
    qtof6600.px.id, nrow(qtof6600.px.fot)
  ),
  Exptype = rep(c('R1', 'R2', 'R3', 'Union'), 3),
  Instrument = rep(c('Exploris480', 'timsTOF', 'QTOF6600'), each=4)
)
barplot_df$Instrument = factor(barplot_df$Instrument, levels = c('QTOF6600', 'timsTOF', 'Exploris480'))
barplot_df$Class = paste(barplot_df$Instrument, barplot_df$Exptype, sep = '_')
barplot_df$Class = factor(barplot_df$Class, levels = barplot_df$Class)

df = barplot_df
p <- ggbarplot(
  df, x = "Class", y = "ID",
  fill = "Instrument", 
  # color = "Exptype", 
  color = "white",
  palette = "jco",
  x.text.angle = 90
)
p + xlab('') + ylab('Protein identification')

jco_col = pal_jco("default", alpha = 1)(10)[c(3,1,2)]
ggdotchart(
  df, x = "Class", y = "ID",
  color = "Instrument",                                # Color by groups
  palette = jco_col, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
  # sorting = "ascending",                        # Sort value in descending order
  add = "segments",                             # Add segments from y = 0 to dots
  ggtheme = theme_pubr(),                        # ggplot2 theme
  group = "Instrument",                                # Order by groups
  dot.size = 6,                                 # Large dot size
  label = round(df$ID),
  add.params = list(color = rep(jco_col, each = 4), size = 2), # Change segment color and size,
  xlab = '',
  ylab = 'Protein identification',
  ylim = c(0, 6000),
  font.label = list(color = "black", size = 10, vjust = -1)
)

# venn
library(VennDiagram)
library(tidyverse)
alpha_per = 0.3
jco_col = pal_jco("default", alpha = 0.5)(10)[c(1)]
venn.diagram(
  x = list(
    R1 = e480.px.col_gps[[1]],
    R2 = e480.px.col_gps[[2]],
    R3 = e480.px.col_gps[[3]]
  ),
  filename = 'temp_figures/03/Manufacturer/03_venn_e480.png',
  main = 'Exploris480',
  # main.just = c(0, 5),
  main.pos = c(0.5,1.05),
  main.cex = 0.5,
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col = jco_col,
  fill = jco_col,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  # cat.pos = c(-2, 27, -135, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
  # cat.col = p_colors
  # rotation = 1,
  # print.mode = 'percent'
)

jco_col = pal_jco("default", alpha = 0.5)(10)[c(2)]
venn.diagram(
  x = list(
    R1 = timstof.px.col_gps[[1]],
    R2 = timstof.px.col_gps[[2]],
    R3 = timstof.px.col_gps[[3]]
  ),
  filename = 'temp_figures/03/Manufacturer/03_venn_timstof.png',
  main = 'timsTOF',
  # main.just = c(0, 5),
  main.pos = c(0.5,1.05),
  main.cex = 0.5,
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col = jco_col,
  fill = jco_col,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  # cat.pos = c(-2, 27, -135, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
  # cat.col = p_colors
  # rotation = 1,
  # print.mode = 'percent'
)

#######################
#Exploris480-timsTOF-QTOF6600
# jco_col = pal_jco("default", alpha = 0.5)(10)[c(1,2)]

jco_col = pal_jco("default", alpha = 1)(10)[c(1,2,3)]
venn.diagram(
  x = list(
    Exploris480 = as.vector(e480.px.fot$GeneSymbol),
    timsTOF = as.vector(timstof.px.fot$GeneSymbol),
    QTOF6600 = as.vector(qtof6600.px.fot$GeneSymbol)
  ),
  filename = 'temp_figures/03/Manufacturer/03_venn_Exploris480_timsTOF_QTOF6600_per.png',
  main = 'QTOF6600 vs Exploris480 vs timsTOF',
  # main.just = c(0, 5),
  main.pos = c(0.5,1.05),
  main.cex = 0.5,
  output = TRUE ,
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col = jco_col,
  fill = jco_col,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.default.pos = "outer",
  # cat.pos = c(-2, 27, -135, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
  # cat.col = p_colors
  # rotation = 1,
  # print.mode = 'percent'
)





#2 identification frequency
e480.row_id = apply(e480.px.fot[,-1], 1, function(x){length(which(x>0))})
timstof.row_id = apply(timstof.px.fot[,-1], 1, function(x){length(which(x>0))})
qtof6600.row_id = apply(qtof6600.px.fot[,-1], 1, function(x){length(which(x>0))})













#3 intensity distribution
get_expr_label <- function(x){
  # x = e480.row_mean
  x.quantile = as.vector(quantile(x, probs = c(1/3, 2/3)))
  x.l_index = which(x<x.quantile[1])
  x.h_index = which(x>x.quantile[2])
  x.expr_label = rep('Medium', length(x))
  x.expr_label[x.l_index] = 'Low'
  x.expr_label[x.h_index] = 'High'
  x.expr_label = factor(x.expr_label, levels = c('Low', 'Medium', 'High'))
  return(x.expr_label)
}

e480.row_mean = rowMeans(e480.px.fot[,-1])
e480.row_expr_label = get_expr_label(e480.row_mean)

timstof.row_mean = rowMeans(timstof.px.fot[,-1])
timstof.row_expr_label = get_expr_label(timstof.row_mean)

qtof6600.row_mean = rowMeans(qtof6600.px.fot[,-1])
qtof6600.row_expr_label = get_expr_label(qtof6600.row_mean)

e480.df = data.frame(
  GeneSymbol = e480.px.fot$GeneSymbol,
  RowID = e480.row_id,
  ExprLab = e480.row_expr_label
)

timstof.df = data.frame(
  GeneSymbol = timstof.px.fot$GeneSymbol,
  RowID = timstof.row_id,
  ExprLab = timstof.row_expr_label
)

qtof6600.df = data.frame(
  GeneSymbol = qtof6600.px.fot$GeneSymbol,
  RowID = qtof6600.row_id,
  ExprLab = qtof6600.row_expr_label
)


table(e480.df$RowID, e480.df$ExprLab)
table(timstof.df$RowID, timstof.df$ExprLab)
table(qtof6600.df$RowID, qtof6600.df$ExprLab)

library(ggforce)
library(ggplot2)
library(webr)
library(moonBook)
df = e480.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
PieDonut_zdd( # 5*5
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(3,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = paste('Exploris480 (Total=', nrow(df), ')', sep = '')
)


df = qtof6600.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
library(ggforce)
PieDonut_zdd(
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(1,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = paste('QTOF6600 (Total=', nrow(df), ')', sep = '')
)

df = timstof.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
library(ggforce)
PieDonut_zdd(
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(3,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = paste('timsTOF (Total=', nrow(df), ')', sep = '')
)





#2 correlation
#Exploris480
if(T){
  df.ref = e480.df
  df = e480.px.fot
  
  cor_method = 'pearson' # spearman
  df.row_id = apply(df[,-1], 1, function(x){length(which(x>0))})
  df = df[df.row_id>=2,]
  df.ref = df.ref[df.row_id>=2,]
  
  px.corv = as.vector(cor(df[,-1], method = cor_method))
  px.corv = unique(px.corv[px.corv<1])
  
  df.l = df[which(df.ref$ExprLab=='Low'),]
  # df.l.row_id = apply(df.l[,-1], 1, function(x){length(which(x>0))})
  # df.l = df.l[df.l.row_id>2,]
  px.corv.l = as.vector(cor(df.l[,-1], method = cor_method))
  px.corv.l = unique(px.corv.l[px.corv.l<1])
  
  df.m = df[which(df.ref$ExprLab=='Medium'),]
  # df.m.row_id = apply(df.m[,-1], 1, function(x){length(which(x>0))})
  # df.m = df.m[df.m.row_id>2,]
  px.corv.m = as.vector(cor(df.m[,-1], method = cor_method))
  px.corv.m = unique(px.corv.m[px.corv.m<1])
  
  df.h = df[which(df.ref$ExprLab=='High'),]
  # df.h.row_id = apply(df.h[,-1], 1, function(x){length(which(x>0))})
  # df.h = df.h[df.h.row_id>2,]
  # df.h.delele.index = which(df.h$GeneSymbol ==	'ACTA1') # contamination
  # df.h = df.h[-df.h.delele.index,]
  px.corv.h = as.vector(cor(df.h[,-1], method = cor_method))
  px.corv.h = unique(px.corv.h[px.corv.h<1])
  
  
  cor_df = data.frame(
    Global = px.corv,
    High = px.corv.h,
    Medium = px.corv.m,
    Low = px.corv.l
  )
  
  df_merge = merge(df.l, df.m, by='GeneSymbol', all = T)
  df_merge = merge(df_merge, df.h, by='GeneSymbol', all = T)
  colnames(df_merge) = c(
    'GeneSymbol',
    'R1_L', 'R2_L', "R3_L",
    'R1_M', 'R2_M', "R3_M",
    'R1_H', 'R2_H', "R3_H"
  )
  df_merge.mat = as.matrix(df_merge[,-1])
  cor(df_merge.mat, use = 'pairwise.complete.obs')
  library('GGally')
  pdf(
    normalizePath(file.path(
      BASE_DIR, 'temp_figures', '03', 'NPB_e480_timstof',
      '03_pairwise_cor_e480.pdf'
    ), mustWork = F), 
    height = 5, width = 5
  )
  # Global
  ggpairs(
    df[,-1],
    columnLabels = c('R1_G', 'R2_G', "R3_G")
  )
  
  # Low
  ggpairs(
    df.l[,-1],
    columnLabels = c('R1_L', 'R2_L', "R3_L")
  )
  
  # Medium
  ggpairs(
    df.m[,-1],
    # lower=list(continuous=my_fn),
    columnLabels = c('R1_M', 'R2_M', "R3_M")
  )
  
  # High
  ggpairs(
    df.h[,-1],
    # lower=list(continuous=my_fn),
    columnLabels = c('R1_H', 'R2_H', "R3_H")
    
  )
  dev.off()
}



if(T){
  #timsTOF
  df.ref = timstof.df
  df = timstof.px.fot
  
  cor_method = 'pearson' # spearman
  df.row_id = apply(df[,-1], 1, function(x){length(which(x>0))})
  df = df[df.row_id>=2,]
  df.ref = df.ref[df.row_id>=2,]
  px.corv = as.vector(cor(df[,-1], method = cor_method))
  px.corv = unique(px.corv[px.corv<1])
  
  df.l = df[which(df.ref$ExprLab=='Low'),]
  px.corv.l = as.vector(cor(df.l[,-1], method = cor_method))
  px.corv.l = unique(px.corv.l[px.corv.l<1])
  
  df.m = df[which(df.ref$ExprLab=='Medium'),]
  px.corv.m = as.vector(cor(df.m[,-1], method = cor_method))
  px.corv.m = unique(px.corv.m[px.corv.m<1])
  
  df.h = df[which(df.ref$ExprLab=='High'),]
  px.corv.h = as.vector(cor(df.h[,-1], method = cor_method))
  px.corv.h = unique(px.corv.h[px.corv.h<1])
  pdf(
    normalizePath(file.path(
      BASE_DIR, 'temp_figures', '03', 'NPB_e480_timstof',
      '03_pairwise_cor_timstof.pdf'
    ), mustWork = F), 
    height = 5, width = 5
  )
  # Global
  ggpairs(
    df[,-1],
    columnLabels = c('R1_G', 'R2_G', "R3_G")
  )
  
  # Low
  ggpairs(
    df.l[,-1],
    columnLabels = c('R1_L', 'R2_L', "R3_L")
  )
  
  # Medium
  ggpairs(
    df.m[,-1],
    # lower=list(continuous=my_fn),
    columnLabels = c('R1_M', 'R2_M', "R3_M")
  )
  
  # High
  ggpairs(
    df.h[,-1],
    # log10(df.h[,-1]+1),
    # lower=list(continuous=my_fn),
    columnLabels = c('R1_H', 'R2_H', "R3_H")
    
  )
  dev.off()
}

if(T){
  #timsTOF
  df.ref = qtof6600.df
  df = qtof6600.px.fot
  
  cor_method = 'pearson' # spearman
  df.row_id = apply(df[,-1], 1, function(x){length(which(x>0))})
  df = df[df.row_id>=2,]
  df.ref = df.ref[df.row_id>=2,]
  px.corv = as.vector(cor(df[,-1], method = cor_method))
  px.corv = unique(px.corv[px.corv<1])
  
  df.l = df[which(df.ref$ExprLab=='Low'),]
  px.corv.l = as.vector(cor(df.l[,-1], method = cor_method))
  px.corv.l = unique(px.corv.l[px.corv.l<1])
  
  df.m = df[which(df.ref$ExprLab=='Medium'),]
  px.corv.m = as.vector(cor(df.m[,-1], method = cor_method))
  px.corv.m = unique(px.corv.m[px.corv.m<1])
  
  df.h = df[which(df.ref$ExprLab=='High'),]
  px.corv.h = as.vector(cor(df.h[,-1], method = cor_method))
  px.corv.h = unique(px.corv.h[px.corv.h<1])
  pdf(
    normalizePath(file.path(
      BASE_DIR, 'temp_figures', '03', 'Manufacturer',
      '03_pairwise_cor_QTOF6600.pdf'
    ), mustWork = F), 
    height = 5, width = 5
  )
  # Global
  ggpairs(
    df[,-1],
    columnLabels = c('R1_G', 'R2_G', "R3_G")
  )
  
  # Low
  ggpairs(
    df.l[,-1],
    columnLabels = c('R1_L', 'R2_L', "R3_L")
  )
  
  # Medium
  ggpairs(
    df.m[,-1],
    # lower=list(continuous=my_fn),
    columnLabels = c('R1_M', 'R2_M', "R3_M")
  )
  
  # High
  ggpairs(
    df.h[,-1],
    # log10(df.h[,-1]+1),
    # lower=list(continuous=my_fn),
    columnLabels = c('R1_H', 'R2_H', "R3_H")
    
  )
  dev.off()
}



#### reproducibility #### 
if(T){
  # reproducibility
  get_px_reproducibility <- function(px_lst, ref_df = NULL, ref_lab = NULL){
    # px_lst = npb.px.col_gps
    # ref_df = npb.df
    # ref_lab = 'Low'
    px_lst.len = length(px_lst)
    ij_repro_list = NULL
    for (i in 1:px_lst.len) {
      i.gps = px_lst[[i]]
      if(!is.null(ref_df) & !is.null(ref_lab)){
        ref_lab.gps = as.vector(ref_df$GeneSymbol[ref_df$ExprLab==ref_lab])
        ref_lab.gps.index = match(ref_lab.gps, i.gps)
        i.gps = i.gps[na.omit(ref_lab.gps.index)]
      }
      
      for (j in 1:px_lst.len) {
        if(i!=j){
          j.gps = px_lst[[j]]
          if(!is.null(ref_df) & !is.null(ref_lab)){
            ref_lab.gps = as.vector(ref_df$GeneSymbol[ref_df$ExprLab==ref_lab])
            ref_lab.gps.index = match(ref_lab.gps, j.gps)
            j.gps = j.gps[na.omit(ref_lab.gps.index)]
          }
          ij_repro = round(
            length(intersect(i.gps, j.gps))/length(j.gps), 2
          )
          ij_repro_list = c(ij_repro_list, ij_repro)
        }
        
      }
    }
    return(ij_repro_list)
    
  } 

  # qtof6600.px.col_gps, qtof6600.df
  qtof6600.repro = get_px_reproducibility(qtof6600.px.col_gps, ref_df = NULL, ref_lab = NULL)
  qtof6600.repro.h = get_px_reproducibility(qtof6600.px.col_gps, ref_df = qtof6600.df, ref_lab = 'High')
  qtof6600.repro.m = get_px_reproducibility(qtof6600.px.col_gps, ref_df = qtof6600.df, ref_lab = 'Medium')
  qtof6600.repro.l = get_px_reproducibility(qtof6600.px.col_gps, ref_df = qtof6600.df, ref_lab = 'Low')
  df.qtof6600 = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(qtof6600.repro.l, qtof6600.repro.m, qtof6600.repro.h, qtof6600.repro)
  )
  df.qtof6600$Intensity = factor(df.qtof6600$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  

  # timstof.px.col_gps, timstof.df
  timstof.repro = get_px_reproducibility(timstof.px.col_gps, ref_df = NULL, ref_lab = NULL)
  timstof.repro.h = get_px_reproducibility(timstof.px.col_gps, ref_df = timstof.df, ref_lab = 'High')
  timstof.repro.m = get_px_reproducibility(timstof.px.col_gps, ref_df = timstof.df, ref_lab = 'Medium')
  timstof.repro.l = get_px_reproducibility(timstof.px.col_gps, ref_df = timstof.df, ref_lab = 'Low')
  df.timstof = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(timstof.repro.l, timstof.repro.m, timstof.repro.h, timstof.repro)
  )
  df.timstof$Intensity = factor(df.timstof$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  
  
  # e480.px.col_gps, e480.df
  e480.repro = get_px_reproducibility(e480.px.col_gps, ref_df = NULL, ref_lab = NULL)
  e480.repro.h = get_px_reproducibility(e480.px.col_gps, ref_df = e480.df, ref_lab = 'High')
  e480.repro.m = get_px_reproducibility(e480.px.col_gps, ref_df = e480.df, ref_lab = 'Medium')
  e480.repro.l = get_px_reproducibility(e480.px.col_gps, ref_df = e480.df, ref_lab = 'Low')
  df.e480 = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(e480.repro.l, e480.repro.m, e480.repro.h, e480.repro)
  )
  df.e480$Intensity = factor(df.e480$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  
  df = rbind.data.frame(
    df.qtof6600, 
    df.timstof,
    df.e480
  )
  df$Instrument = rep(
    c('QTOF6600', 'timsTOF', 'Exploris480'), each = 24
  )
  df$Instrument = factor(df$Instrument, levels = c('QTOF6600', 'timsTOF', 'Exploris480'))
  
  # '#33C860', # 绿色
  # '#81B0FF', # 蓝色
  # '#F9918A'  # 红色
  p <- ggboxplot(
    df, x = "Intensity", y = "Reproducibity",
    color = "Intensity", palette =c("#33C860", "#81B0FF", "#F9918A", '#EFC000FF'),
    add = "jitter", shape = "Instrument", facet.by = 'Instrument'
  )
  p
  
  
  
  
  
  
  # merge cor
  # qtof6600.px.fot
  # e480.px.fot
  if(T){
    qtof6600.px.fot.row_id = apply(qtof6600.px.fot[,-1], 1, function(x){length(which(x>0))})
    timstof.px.fot.row_id = apply(timstof.px.fot[,-1], 1, function(x){length(which(x>0))})
    e480.px.fot.row_id = apply(e480.px.fot[,-1], 1, function(x){length(which(x>0))})
    
    px.fot = merge(
      qtof6600.px.fot[qtof6600.px.fot.row_id>=1,],
      merge(
        timstof.px.fot[timstof.px.fot.row_id>=1,], 
        e480.px.fot[e480.px.fot.row_id>=1,], 
        by = 'GeneSymbol', all = T
      ),
      by = 'GeneSymbol', all = T
    )
    px.fot.mat = as.matrix(px.fot[,-1])
    colnames(px.fot.mat) = c(
      paste('QTOF6600', c('R1', 'R2', 'R3'), sep = '_'),
      paste('timsTOF', c('R1', 'R2', 'R3'), sep = '_'),
      paste('Exploris480', c('R1', 'R2', 'R3'), sep = '_')
    )
    na.index = which(is.na(px.fot.mat))
    px.fot.mat[na.index] = 0
    cor_method = 'pearson'
    cormat <- round(cor(px.fot.mat, method = cor_method),2)
    library(reshape2)
    melted_cormat <- melt(cormat)
    library(ggplot2) # 6*4.5
    p1 <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) 
    p1 <- p1 +  geom_tile() 
    p1 <- p1 + xlab('Site-Replicates') + ylab('Site-Replicates')
    p1 <- p1 + theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
    # p1 <- p1 + scale_fill_gradient2(
    #   low = "#132B43", high ="#56B1F7", mid = "white",
    #   limits = c(0.5, 1), breaks = seq(0.5,1,0.1),
    #   guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)
    # )
    p1 <- p1 + scale_fill_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)
    ) + labs(fill = "PCorr")
    p1
    
    # col<- colorRampPalette(c("blue", "white", "red"))(20)
    # heatmap(x = cormat, col = col, symm = TRUE)
  }
}







# High
e480.h.gps = e480.df$GeneSymbol[e480.df$ExprLab=='High']
timstof.h.gps = timstof.df$GeneSymbol[timstof.df$ExprLab=='High']
length(intersect(e480.h.gps, timstof.h.gps))/length(union(e480.h.gps, timstof.h.gps))

# Medium
e480.m.gps = e480.df$GeneSymbol[e480.df$ExprLab=='Medium']
timstof.m.gps = timstof.df$GeneSymbol[timstof.df$ExprLab=='Medium']
length(intersect(e480.m.gps, timstof.m.gps))/length(union(e480.m.gps, timstof.m.gps))

# Low
e480.l.gps = e480.df$GeneSymbol[e480.df$ExprLab=='Low']
timstof.l.gps = timstof.df$GeneSymbol[timstof.df$ExprLab=='Low']
length(intersect(e480.l.gps, timstof.l.gps))/length(union(e480.l.gps, timstof.l.gps))


# Sankey
e480_timstof.sankey_df <- generate_sankey_df(e480.df, timstof.df, c('Exploris480', 'timsTOF'))
links_df.tmp = e480_timstof.sankey_df

links_df = data.frame(
  source = paste(links_df.tmp$MS1, links_df.tmp$Lab1, sep='-'),
  target = paste(links_df.tmp$MS2, links_df.tmp$Lab2, sep='-'),
  value = links_df.tmp$Value
)
links_df$group = links_df$target

links_df$source = factor(
  links_df$source, levels = c('Exploris480-High', 'Exploris480-Medium', 'Exploris480-Low', 'Exploris480-NA')
)
links_df$target = factor(
  links_df$target, levels = c('timsTOF-High', 'timsTOF-Medium', 'timsTOF-Low', 'timsTOF-NA')
)

links = links_df[order(links_df$source),]
# links = links[c(1,3,2,4, 5,7,6,8, 9,11,10,12, 13,15,14),]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
# '#33C860', # 绿色
# '#81B0FF', # 蓝色
# '#F9918A'  # 红色
my_color <- 'd3.scaleOrdinal() .domain(["Exploris480-High", "Exploris480-Medium", "Exploris480-Low", "Exploris480-NA", "timsTOF-High", "timsTOF-Medium", "timsTOF-Low", "timsTOF-NA"]) .range(["#F9918A", "#81B0FF", "#33C860", "grey", "#F9918A", "#81B0FF", "#33C860", "grey"])'
my_color <- 'd3.scaleOrdinal() .domain(["Exploris480-High", "Exploris480-Medium", "Exploris480-Low", "Exploris480-NA", "timsTOF-High", "timsTOF-Medium", "timsTOF-Low", "timsTOF-NA"]) .range(["#F9918A", "#81B0FF", "#33C860", "grey", "#F9918A", "#81B0FF", "#33C860", "grey"])'
my_color <- 'd3.scaleOrdinal() .domain(["Exploris480-High", "Exploris480-Medium", "Exploris480-Low", "Exploris480-NA", "timsTOF-High", "timsTOF-Medium", "timsTOF-Low", "timsTOF-NA"]) .range(["#F9918A", "#81B0FF", "#33C860", "grey", "#F9918A", "#81B0FF", "#33C860", "grey"])'
p <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", 
  NodeID = "name", 
  fontSize = 15,
  colourScale = my_color,
  LinkGroup="group",
  iterations = 0,
  sinksRight = FALSE
)
p




# Sankey
qtof6600_e480.sankey_df <- generate_sankey_df(qtof6600.df, e480.df, c('QTOF6600', 'Exploris480'))
links_df.tmp = qtof6600_e480.sankey_df

links_df = data.frame(
  source = paste(links_df.tmp$MS1, links_df.tmp$Lab1, sep='-'),
  target = paste(links_df.tmp$MS2, links_df.tmp$Lab2, sep='-'),
  value = links_df.tmp$Value
)
links_df$group = links_df$target

links_df$source = factor(
  links_df$source, levels = c('QTOF6600-High', 'QTOF6600-Medium', 'QTOF6600-Low', 'QTOF6600-NA')
)
links_df$target = factor(
  links_df$target, levels = c('Exploris480-High', 'Exploris480-Medium', 'Exploris480-Low', 'Exploris480-NA')
)

links = links_df[order(links_df$source),]
# links = links[c(1,3,2,4, 5,7,6,8, 9,11,10,12, 13,15,14),]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
# '#33C860', # 绿色
# '#81B0FF', # 蓝色
# '#F9918A'  # 红色
my_color <- 'd3.scaleOrdinal() .domain(["QTOF6600-High", "QTOF6600-Medium", "QTOF6600-Low", "QTOF6600-NA", "Exploris480-High", "Exploris480-Medium", "Exploris480-Low", "Exploris480-NA"]) .range(["#F9918A", "#81B0FF", "#33C860", "grey", "#F9918A", "#81B0FF", "#33C860", "grey"])'
p <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", 
  NodeID = "name", 
  fontSize = 15,
  colourScale = my_color,
  LinkGroup="group",
  iterations = 0,
  sinksRight = FALSE
)
p



# Sankey
qtof6600_timsTOF.sankey_df <- generate_sankey_df(qtof6600.df, timstof.df, c('QTOF6600', 'timsTOF'))
links_df.tmp = qtof6600_timsTOF.sankey_df

links_df = data.frame(
  source = paste(links_df.tmp$MS1, links_df.tmp$Lab1, sep='-'),
  target = paste(links_df.tmp$MS2, links_df.tmp$Lab2, sep='-'),
  value = links_df.tmp$Value
)
links_df$group = links_df$target

links_df$source = factor(
  links_df$source, levels = c('QTOF6600-High', 'QTOF6600-Medium', 'QTOF6600-Low', 'QTOF6600-NA')
)
links_df$target = factor(
  links_df$target, levels = c('timsTOF-High', 'timsTOF-Medium', 'timsTOF-Low', 'timsTOF-NA')
)

links = links_df[order(links_df$source),]
# links = links[c(1,3,2,4, 5,7,6,8, 9,11,10,12, 13,15,14),]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
# '#33C860', # 绿色
# '#81B0FF', # 蓝色
# '#F9918A'  # 红色
my_color <- 'd3.scaleOrdinal() .domain(["QTOF6600-High", "QTOF6600-Medium", "QTOF6600-Low", "QTOF6600-NA", "timsTOF-High", "timsTOF-Medium", "timsTOF-Low", "timsTOF-NA"]) .range(["#F9918A", "#81B0FF", "#33C860", "grey", "#F9918A", "#81B0FF", "#33C860", "grey"])'
p <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", 
  NodeID = "name", 
  fontSize = 15,
  colourScale = my_color,
  LinkGroup="group",
  iterations = 0,
  sinksRight = FALSE
)
p































# 相关性
qtof6600.df$FOT = apply(qtof6600.px.fot[,-1], 1, mean)
timstof.df$FOT = apply(timstof.px.fot[,-1], 1, mean)
e480.df$FOT = apply(e480.px.fot[,-1], 1, mean)



merge_df = merge(
  qtof6600.df,
  merge(timstof.df, e480.df, by = 'GeneSymbol', all = T),
  by = 'GeneSymbol', all = T
)

merge_df = na.omit(merge_df)

# Global
cor(merge_df$FOT.x, merge_df$FOT.y)
cor.g <- ggpairs(
  log10(merge_df[,c(4,7)]),
  columnLabels = c('Exploris480-Global', 'timsTOF-Global'),
  lower=list(continuous=my_fn)
)


# High
merge_df.high = merge_df[which(
  merge_df$ExprLab.x=='High' & merge_df$ExprLab.y=='High' & merge_df$ExprLab =='High'
),]
cor(merge_df.high$FOT.x, merge_df.high$FOT.y)
cor.h <- ggpairs(
  log10(merge_df.high[,c(4,7)]),
  columnLabels = c('Exploris480-High', 'timsTOF-High')
)

# Medium
merge_df.medium = merge_df[which(
  merge_df$ExprLab.x=='Medium' & merge_df$ExprLab.y=='Medium'
),]
cor(merge_df.medium$FOT.x, merge_df.medium$FOT.y)
cor.m <- ggpairs(
  log10( merge_df.medium[,c(4,7)]),
  columnLabels = c('Exploris480-Medium', 'timsTOF-Medium')
)

# Low
merge_df.low = merge_df[which(
  merge_df$ExprLab.x=='Low' & merge_df$ExprLab.y=='Low'
),]
cor(merge_df.low$FOT.x, merge_df.low$FOT.y)
cor.l <- ggpairs(
  log10( merge_df.low[,c(4,7)]),
  columnLabels = c('Exploris480-Low', 'timsTOF-Low')
)


get_ggplot_cor_density <- function(dat, label){
  # dat = log10(dat)
  df <- data.frame(
    x = dat[,1],
    y = dat[,2]
  )
  df.cor = cor(dat$FOT.x, dat$FOT.y)
  df.cor = round(df.cor, 2)
  df.cor.text = paste(label, '(N=', nrow(df), ')', ': Cor=', df.cor, sep = '')
  
  df = log10(df)
  p <- ggplot(df, aes(x=x, y=y))
  p <- p + stat_density2d(aes(fill=..density..), geom="tile", contour = F)
  p <- p + scale_fill_gradientn(colours=magma(100)[10:100])
  p <- p + xlab("log10(FOT)/Exploris480") + ylab("log10(FOT)/timsTOF") 
  p <- p + scale_x_continuous(breaks=seq(-4, 4, 1)) + scale_y_continuous(breaks=seq(-4, 4, 1))
  p <- p + xlim(-4, 4) + ylim(-4,4)
  p <- p + theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )
  p <- p + annotate("text", x = 0, y = -3.5, label = df.cor.text, color='yellow', size=3)
  return(p)
}

p.g = get_ggplot_cor_density((merge_df[,c(4,7)]), 'Global')
p.h = get_ggplot_cor_density((merge_df.high[,c(4,7)]), 'High')
p.m = get_ggplot_cor_density((merge_df.medium[,c(4,7)]), 'Medium')
p.l = get_ggplot_cor_density((merge_df.low[,c(4,7)]), 'Low')

ggarrange(#8*6
  p.l, p.m, p.h, p.g, 
  nrow = 2, ncol = 2,
  labels = c('A', 'B', 'C', 'D')
)


# CV
common_genes = as.vector(merge_df$GeneSymbol)
h_index = which(merge_df$ExprLab.x=='High' & merge_df$ExprLab.y=='High' & merge_df$ExprLab=='High')
m_index = which(merge_df$ExprLab.x=='Medium' & merge_df$ExprLab.y=='Medium' & merge_df$ExprLab=='Medium')
l_index = which(merge_df$ExprLab.x=='Low' & merge_df$ExprLab.y=='Low' & merge_df$ExprLab=='Low')

qtof6600.px.fot.common = qtof6600.px.fot[match(common_genes,qtof6600.px.fot$GeneSymbol),]
timstof.px.fot.common = timstof.px.fot[match(common_genes,timstof.px.fot$GeneSymbol),]
e480.px.fot.common = e480.px.fot[match(common_genes,e480.px.fot$GeneSymbol),]


get_cv_list <- function(df){
  # df = e480.px.fot.common
  df.row_cv = apply(df[,-1], 1, function(x){
    sd(x)/mean(x)
  })
  df.row_cv = round(df.row_cv, 2)
  return(df.row_cv)
}

qtof6600.g.cv = get_cv_list(qtof6600.px.fot.common)
qtof6600.h.cv = get_cv_list(qtof6600.px.fot.common[merge_df$ExprLab=='High',])
qtof6600.m.cv = get_cv_list(qtof6600.px.fot.common[merge_df$ExprLab=='Medium',])
qtof6600.l.cv = get_cv_list(qtof6600.px.fot.common[merge_df$ExprLab=='Low',])


timstof.g.cv = get_cv_list(timstof.px.fot.common)
timstof.h.cv = get_cv_list(timstof.px.fot.common[merge_df$ExprLab.x=='High',])
timstof.m.cv = get_cv_list(timstof.px.fot.common[merge_df$ExprLab.x=='Medium',])
timstof.l.cv = get_cv_list(timstof.px.fot.common[merge_df$ExprLab.x=='Low',])

e480.g.cv = get_cv_list(e480.px.fot.common)
e480.h.cv = get_cv_list(e480.px.fot.common[merge_df$ExprLab.y=='High',])
e480.m.cv = get_cv_list(e480.px.fot.common[merge_df$ExprLab.y=='Medium',])
e480.l.cv = get_cv_list(e480.px.fot.common[merge_df$ExprLab.y=='Low',])



timstofX_timstof.g.cv = get_cv_list(
  data.frame(e480.px.fot.common, timstof.px.fot.common[,-1])
)
timstofX_timstof.h.cv = get_cv_list(
  data.frame(e480.px.fot.common, timstof.px.fot.common[,-1])[h_index,]
)
timstofX_timstof.m.cv = get_cv_list(
  data.frame(e480.px.fot.common, timstof.px.fot.common[,-1])[m_index,]
)
timstofX_timstof.l.cv = get_cv_list(
  data.frame(e480.px.fot.common, timstof.px.fot.common[,-1])[l_index,]
)

bp <- boxplot(
  qtof6600.l.cv, qtof6600.m.cv, qtof6600.h.cv, qtof6600.g.cv,
  timstof.l.cv, timstof.m.cv, timstof.h.cv, timstof.g.cv,
  e480.l.cv, e480.m.cv, e480.h.cv, e480.g.cv
)

# '#33C860', # 绿色
# '#81B0FF', # 蓝色
# '#F9918A'  # 红色

if(T){#8*6
  bp <- boxplot(
    NA, 
    xlim = c(0, 12.5), ylim = c(0, 3), ylab = 'Coefficient of Variation',
    xaxt = 'n', yaxt = 'n'
  )
  
  boxplot(
    qtof6600.l.cv, qtof6600.m.cv, qtof6600.h.cv, qtof6600.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25), boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  boxplot(
    timstof.l.cv, timstof.m.cv, timstof.h.cv, timstof.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25)+4.5, boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  boxplot(
    e480.l.cv, e480.m.cv, e480.h.cv, e480.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25)+4.5*2, boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  abline(v = 4, lty = 3, lwd = 2, col = 'grey')
  abline(v = 8.5, lty = 3, lwd = 2, col = 'grey')
  axis(
    1, at = c(1.75, 6.25, 10.75), labels = c('QTOF6600', 'timsTOF', 'Exploris480')
  )
  
  axis(
    2, at = seq(0, 3, 0.5), 
    labels = c('0', '50%', '100%', '150%', '200%', '250%', '300%'), 
    las = 1, cex.axis = 0.8
  )
  axis(
    2, at = c(0.2, 0.3), 
    labels = c('20%', '30%'), 
    tcl = -0.3, cex.axis = 0.6, las = 1, col.axis = 'red', col = 'red'
  )
  abline(h = c(0.2, 0.3), lty = 3, lwd = 2, col = 'red')
}

# DEPs
# meta_data.NPB.e480 
# meta_data.NPB.timstof

# D5
e480.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.e480$Exp_code)
timstof.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.timstof$Exp_code)

get_deps <- function(df, group, min_fc=2){
  # df = e480.fot
  # group = meta_data.NPB.e480$UPerson
  
  df.genes = as.vector(unlist(df[,1]))
  df.mat = as.matrix(df[,-1])
  group.uid = unique(group)
  group.uid.count = length(group.uid)
  i.end = group.uid.count - 1
  i_j_id = NULL
  deps_list = list()
  deps_list.index = 0
  
  for (i in 1:i.end) {
    ctrl.uid = group.uid[i]
    ctrl.index = which(group==ctrl.uid)
    j.start = i + 1
    for (j in j.start:group.uid.count) {
      case.uid = group.uid[j]
      case.index = which(group==case.uid) 
      i_j_id = c(i_j_id, paste(ctrl.uid, case.uid, sep = ' VS '))
      
      ctrl.mean = rowMeans(df.mat[,ctrl.index])
      case.mean = rowMeans(df.mat[,case.index])
      # pval = apply(data.frame(df.mat[,c(ctrl.index, case.index)]), 1, function(x){
      #   wilcox.test(x[1:3], x[4:6])$p.value
      # })
      fc = ctrl.mean/case.mean
      fc.index = which(
        !is.nan(fc) & (fc<1/min_fc | fc>min_fc)#  & pval<0.05
      )
      i_j_df = data.frame(
        GeneSymbol = df.genes,
        FC = fc,
        Ctrl_id = apply(df.mat[,ctrl.index], 1, function(x){length(which(x>0))}),
        Case_id = apply(df.mat[,case.index], 1, function(x){length(which(x>0))})
      )[fc.index,]
      deps_list.index = deps_list.index + 1
      deps_list[[deps_list.index]] = i_j_df
    }
  }
  names(deps_list) = i_j_id
  return(deps_list)
}

intersect_genes = intersect(as.vector(e480.fot$GeneSymbol), as.vector(timstof.fot$GeneSymbol))
min_fc = 2
e480.fot.deps_list = get_deps(e480.fot[match(intersect_genes, e480.fot$GeneSymbol),], meta_data.NPB.e480$UPerson, min_fc=2)
timstof.fot.deps_list = get_deps(timstof.fot[match(intersect_genes, timstof.fot$GeneSymbol),], meta_data.NPB.timstof$UPerson, min_fc=2)

get_consensus_deps = function(df1, df2, min_fc = 2){
  df1 = e480.fot.deps_list[[1]]
  df2 = timstof.fot.deps_list[[1]]
  
  # df1 = df1[which(
  #   df1$Ctrl_id>1 | df1$Case_id>1
  # ),]
  print(dim(df1))
  
  # df2 = df2[which(
  #   df2$Ctrl_id>1 | df2$Case_id>1
  # ),]
  print(dim(df2))
  
  c_gps.h = intersect(
    df1$GeneSymbol[df1$FC>min_fc], 
    df2$GeneSymbol[df2$FC>min_fc]
  )
  c_gps.l = intersect(
    df1$GeneSymbol[df1$FC<min_fc], 
    df2$GeneSymbol[df2$FC<1/min_fc]
  )
  
  c_gps = c(c_gps.h, c_gps.l)
  
  print(length(c_gps.h))
  print(length(c_gps.l))
  
  
  
  df1$Label = rep('Non', nrow(df1))
  df1$Label[df1$FC>min_fc] = 'UP'
  df1$Label[df1$FC<1/min_fc] = 'DOWN'
  
  df2$Label = rep('Non', nrow(df2))
  df2$Label[df2$FC>min_fc] = 'UP'
  df2$Label[df2$FC<1/min_fc] = 'DOWN'
  
}


if(T){
  # Sankey
  sankey_df = merge(df1, df2, by = 'GeneSymbol', all = T)
  sankey_df = sankey_df[,c(1,5,9)]
  
  
  colnames(sankey_df) = c('GeneSymbol', 'source', 'target')
  st_labs_list = NULL
  for (i in 1:nrow(sankey_df)) {
    s_lab = as.vector(sankey_df[i,2])
    t_lab = as.vector(sankey_df[i,3])
    if(is.na(s_lab)){
      s_lab = 'NA'
    }
    if(is.na(t_lab)){
      t_lab = 'NA'
    }
    s_lab = paste('Exploris480', '-', s_lab, sep = '')
    t_lab = paste('timsTOF', '-', t_lab, sep = '')
    st_labs = paste(s_lab, t_lab, sep = '_')
    st_labs_list = c(st_labs_list, st_labs)
  }
  st_labs_list.table = table(st_labs_list)
  st_labs_list.table_names = names(st_labs_list.table)
  links_df = NULL
  for (i in 1:length(st_labs_list.table_names)) {
    x = st_labs_list.table_names[i]
    x = strsplit(x, split = '_')[[1]]
    links_df = rbind(links_df, x)
  }
  links_df = as.data.frame(links_df)
  rownames(links_df) = st_labs_list.table_names
  colnames(links_df) = c('source', 'target')
  links_df$value = as.vector(st_labs_list.table)
  links_df$group = links_df$target
  
  
  links_df$source = factor(
    links_df$source, levels = c('Exploris480-DOWN', 'Exploris480-NA', 'Exploris480-UP')
  )
  links_df$target = factor(
    links_df$target, levels = c('timsTOF-UP', 'timsTOF-DOWN', 'timsTOF-NA')
  )
  
  links = links_df[order(links_df$source),]
  # links = links[c(3,1,2, 6,4,5, 8,7),]
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )
  # nodes$name = factor(nodes$name, levels = c('Exploris480-UP', 'Exploris480-DOWN', 'Exploris480-NA', 'timsTOF-UP', 'timsTOF-DOWN', 'timsTOF-NA'))
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Make the Network
  # '#33C860', # 绿色
  # '#81B0FF', # 蓝色
  # '#F9918A'  # 红色
  my_color <- 'd3.scaleOrdinal() .domain(["Exploris480-UP", "Exploris480-DOWN", "Exploris480-NA", "timsTOF-UP", "timsTOF-DOWN", "timsTOF-NA"]) .range(["#F9918A",  "#33C860", "grey", "#F9918A",  "#33C860", "grey"])'
  p <- sankeyNetwork(
    Links = links, 
    Nodes = nodes,
    Source = "IDsource", Target = "IDtarget",
    Value = "value", 
    NodeID = "name", 
    fontSize = 15,
    colourScale = my_color,
    LinkGroup="group",
    sinksRight = FALSE
  )
  p
}





if(T){
  # PCA
  qtof6600.fot = get_fot_df_by_expnames(fot_data, meta_data.BRK.QTOF6600$Exp_code)
  timstof.fot = get_fot_df_by_expnames(fot_data, meta_data.BRK.timsTOF$Exp_code)
  e480.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.E480$Exp_code)
  
  # qtof6600.fot
  qtof6600.fot.row_mean = rowMeans(qtof6600.fot[,-1])
  qtof6600.fot.row_expr_label = get_expr_label(qtof6600.fot.row_mean)
  qtof6600.fot.df = data.frame(
    GeneSymbol = qtof6600.fot$GeneSymbol,
    RowID = apply(qtof6600.fot[,-1], 1, function(x){length(which(x>0))}),
    MeanFOT = qtof6600.fot.row_mean,
    ExprLab = qtof6600.fot.row_expr_label
  )
  
  # timstof.fot
  timstof.fot.row_mean = rowMeans(timstof.fot[,-1])
  timstof.fot.row_expr_label = get_expr_label(timstof.fot.row_mean)
  timstof.fot.df = data.frame(
    GeneSymbol = timstof.fot$GeneSymbol,
    RowID = apply(timstof.fot[,-1], 1, function(x){length(which(x>0))}),
    MeanFOT = timstof.fot.row_mean,
    ExprLab = timstof.fot.row_expr_label
  )
  
  # e480.fot
  e480.fot.row_mean = rowMeans(e480.fot[,-1])
  e480.fot.row_expr_label = get_expr_label(e480.fot.row_mean)
  e480.fot.df = data.frame(
    GeneSymbol = e480.fot$GeneSymbol,
    RowID = apply(e480.fot[,-1], 1, function(x){length(which(x>0))}),
    MeanFOT = e480.fot.row_mean,
    ExprLab = e480.fot.row_expr_label
  )
  
  # Single PCA
  # QTOF6600
  if(T){
    fot_df = qtof6600.fot
    fot_df.mat = as.matrix(fot_df[,-1])
    fot_df.mat.scale = t(scale(t(fot_df.mat)))
    fot_df.ref = qtof6600.fot.df
    fot_group = rep(p_flags, each = 3)
    par(mfrow = c(2,2))
    if(T){
      # Low
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('QTOF6600; Low-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main,
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    
    if(T){
      # Medium
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Medium',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('QTOF6600; Medium-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    if(T){
      # High
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='High',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('QTOF6600; High-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    if(T){
      # Global
      mat = fot_df.mat.scale #[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('QTOF6600; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
  }
  
  # timsTOF
  if(T){
    fot_df = timstof.fot
    fot_df.ref = timstof.fot.df
    fot_df.mat = as.matrix(fot_df[,-1])
    fot_df.mat.scale = t(scale(t(fot_df.mat)))
    fot_group = rep(p_flags, each = 3)
    par(mfrow = c(2,2))
    if(T){
      # Low
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('timsTOF; Low-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      sub = 'Exploris480-Exploris480'
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main,
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    
    if(T){
      # Medium
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Medium',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('timsTOF; Medium-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    if(T){
      # High
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='High',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('timsTOF; High-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    if(T){ # 8*8
      # Global
      mat = fot_df.mat.scale #[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('timsTOF; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
  }
  
  
  
  # Exploris480
  if(T){
    fot_df = e480.fot
    fot_df.ref = e480.fot.df
    fot_df.mat = as.matrix(fot_df[,-1])
    fot_df.mat.scale = t(scale(t(fot_df.mat)))
    fot_group = rep(p_flags, each = 3)
    par(mfrow = c(2,2))
    if(T){
      # Low
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('Exploris480; Low-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      sub = 'Exploris480-Exploris480'
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main,
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    
    if(T){
      # Medium
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Medium',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('Exploris480; Medium-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    if(T){
      # High
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='High',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('Exploris480; High-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
    
    
    if(T){ # 8*8
      # Global
      mat = fot_df.mat.scale #[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste('Exploris480; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
  }
}
