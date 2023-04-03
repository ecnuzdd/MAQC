# fot_data
# meta_data

meta_data.NPB.QEPlus = meta_data[which(
  meta_data$Labcode == 'NPB' 
  & meta_data$SInstrument == 'QE-Plus'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]

meta_data.NPB.QEHF = meta_data[which(
  meta_data$Labcode == 'NPB' 
  & meta_data$SInstrument == 'QE-HF'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]


# D5
px_list = c('D5', 'D6', 'F7', 'M8')
px = px_list[1]
qeplus.px.index = which(meta_data.NPB.QEPlus$UPerson == px)
qeplus.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEPlus$Exp_code[qeplus.px.index])
# delele.index = which(qeplus.px.fot$GeneSymbol ==	'ACTA1') # contamination
# if(length(which(delele.index>0))){
#   print(delele.index)
#   qeplus.px.fot = qeplus.px.fot[-delele.index,] # 峰度排名在单个实验排名第一，但只被鉴定到一次
# }

qehf.px.index = which(meta_data.NPB.QEHF$UPerson == px)
qehf.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEHF$Exp_code[qehf.px.index])

#1 identification
qeplus.px.col_gps = apply(qeplus.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(qeplus.px.fot[,1])))
qehf.px.col_gps = apply(qehf.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(qehf.px.fot[,1])))


qeplus.px.id = as.vector(unlist(lapply(qeplus.px.col_gps, length)))
qehf.px.id = as.vector(unlist(lapply(qehf.px.col_gps, length)))
barplot_df = data.frame(
  ID =  c(qeplus.px.id, nrow(qeplus.px.fot), qehf.px.id, nrow(qehf.px.fot)),
  Exptype = rep(c('R1', 'R2', 'R3', 'Union'), 2),
  Instrument = rep(c('QE-Plus', 'QE-HF'), each=4)
)
barplot_df$Instrument = factor(barplot_df$Instrument, levels = c('QE-Plus', 'QE-HF'))
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


ggdotchart(
  df, x = "Class", y = "ID",
  color = "Instrument",                                # Color by groups
  palette = c("#0073C2FF", "#EFC000FF"), # Custom color palette
  # sorting = "ascending",                        # Sort value in descending order
  add = "segments",                             # Add segments from y = 0 to dots
  ggtheme = theme_pubr(),                        # ggplot2 theme
  group = "Instrument",                                # Order by groups
  dot.size = 6,                                 # Large dot size
  label = round(df$ID),
  add.params = list(color = rep(c("#0073C2FF", "#EFC000FF"), each = 4), size = 2), # Change segment color and size,
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
    R1 = qeplus.px.col_gps[[1]],
    R2 = qeplus.px.col_gps[[2]],
    R3 = qeplus.px.col_gps[[3]]
  ),
  filename = 'Venn/03_venn_QE_Plus.png',
  main = 'QE-Plus',
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
    R1 = qehf.px.col_gps[[1]],
    R2 = qehf.px.col_gps[[2]],
    R3 = qehf.px.col_gps[[3]]
  ),
  filename = 'Venn/03_venn_QE_HF.png',
  main = 'QE-HF',
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

#QE-Plus-QE-HF
jco_col = pal_jco("default", alpha = 0.5)(10)[c(1,2)]
venn.diagram(
  x = list(
    QE_Plus = as.vector(qeplus.px.fot$GeneSymbol),
    QE_HF = as.vector(qehf.px.fot$GeneSymbol)
  ),
  filename = 'Venn/03_venn_QE-Plus_QE-HF_per.png',
  main = 'QE-Plus vs QE-HF',
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





#2 identification frequency
QE_Plus.row_id = apply(qeplus.px.fot[,-1], 1, function(x){length(which(x>0))})
QE_HF.row_id = apply(qehf.px.fot[,-1], 1, function(x){length(which(x>0))})













#3 intensity distribution
get_expr_label <- function(x){
  # x = QE_Plus.row_mean
  x.quantile = as.vector(quantile(x, probs = c(1/3, 2/3)))
  x.l_index = which(x<x.quantile[1])
  x.h_index = which(x>x.quantile[2])
  x.expr_label = rep('Medium', length(x))
  x.expr_label[x.l_index] = 'Low'
  x.expr_label[x.h_index] = 'High'
  x.expr_label = factor(x.expr_label, levels = c('Low', 'Medium', 'High'))
  return(x.expr_label)
}

QE_Plus.row_mean = rowMeans(qeplus.px.fot[,-1])
QE_Plus.row_expr_label = get_expr_label(QE_Plus.row_mean)

QE_HF.row_mean = rowMeans(qehf.px.fot[,-1])
QE_HF.row_expr_label = get_expr_label(QE_HF.row_mean)

QE_Plus.df = data.frame(
  GeneSymbol = qeplus.px.fot$GeneSymbol,
  RowID = QE_Plus.row_id,
  ExprLab = QE_Plus.row_expr_label
)

QE_HF.df = data.frame(
  GeneSymbol = qehf.px.fot$GeneSymbol,
  RowID = QE_HF.row_id,
  ExprLab = QE_HF.row_expr_label
)


table(QE_Plus.df$RowID, QE_Plus.df$ExprLab)
table(QE_HF.df$RowID, QE_HF.df$ExprLab)
library(ggforce)
library(ggplot2)
library(webr)
library(moonBook)
df = QE_Plus.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
PieDonut_zdd( # 5*5
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(1,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = paste('QE-Plus (Total=', nrow(df), ')', sep = '')
)


df = QE_HF.df
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
  title = paste('QE-HF (Total=', nrow(df), ')', sep = '')
)






#2 correlation
#QE-Plus
if(T){
  df.ref = QE_Plus.df
  df = qeplus.px.fot
  
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
      BASE_DIR, 'temp_figures', '03', 'NPB_QEPlus_QEHF',
      '03_pairwise_cor_QE_Plus.pdf'
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
  #QE-HF
  df.ref = QE_HF.df
  df = qehf.px.fot
  
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
      BASE_DIR, 'temp_figures', '03', 'NPB_QEPlus_QEHF',
      '03_pairwise_cor_QE_HF.pdf'
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



# High
QE_Plus.h.gps = QE_Plus.df$GeneSymbol[QE_Plus.df$ExprLab=='High']
QE_HF.h.gps = QE_HF.df$GeneSymbol[QE_HF.df$ExprLab=='High']
length(intersect(QE_Plus.h.gps, QE_HF.h.gps))/length(union(QE_Plus.h.gps, QE_HF.h.gps))

# Medium
QE_Plus.m.gps = QE_Plus.df$GeneSymbol[QE_Plus.df$ExprLab=='Medium']
QE_HF.m.gps = QE_HF.df$GeneSymbol[QE_HF.df$ExprLab=='Medium']
length(intersect(QE_Plus.m.gps, QE_HF.m.gps))/length(union(QE_Plus.m.gps, QE_HF.m.gps))

# Low
QE_Plus.l.gps = QE_Plus.df$GeneSymbol[QE_Plus.df$ExprLab=='Low']
QE_HF.l.gps = QE_HF.df$GeneSymbol[QE_HF.df$ExprLab=='Low']
length(intersect(QE_Plus.l.gps, QE_HF.l.gps))/length(union(QE_Plus.l.gps, QE_HF.l.gps))



# Sankey
sankey_df = merge(QE_Plus.df, QE_HF.df, by = 'GeneSymbol', all = T)
sankey_df = sankey_df[,c(1,3,5)]
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
  s_lab = paste('QE-Plus', '-', s_lab, sep = '')
  t_lab = paste('QE-HF', '-', t_lab, sep = '')
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
# links_df$group = factor(
#   links_df$group, levels = c(
#     st_labs_list.table_names[c(1,3,2,4, 9,11,10,12, 5,7,6,8, 13,15,14)]
#   )
# )

links_df$source = factor(
  links_df$source, levels = c('QE-Plus-High', 'QE-Plus-Medium', 'QE-Plus-Low', 'QE-Plus-NA')
)
links_df$target = factor(
  links_df$target, levels = c('QE-HF-High', 'QE-HF-Medium', 'QE-HF-Low', 'QE-HF-NA')
)

links = links_df[order(links_df$source),]
links = links[c(1,3,2,4, 5,7,6,8, 9,11,10,12, 13,15,14),]

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
my_color <- 'd3.scaleOrdinal() .domain(["QE-Plus-High", "QE-Plus-Medium", "QE-Plus-Low", "QE-Plus-NA", "QE-HF-High", "QE-HF-Medium", "QE-HF-Low", "QE-HF-NA"]) .range(["#F9918A", "#81B0FF", "#33C860", "grey", "#F9918A", "#81B0FF", "#33C860", "grey"])'
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
QE_Plus.df$FOT = apply(qeplus.px.fot[,-1], 1, mean)
QE_HF.df$FOT = apply(qehf.px.fot[,-1], 1, mean)
merge_df = merge(QE_Plus.df, QE_HF.df, by = 'GeneSymbol', all = T)
merge_df = na.omit(merge_df)

# Global
cor(merge_df$FOT.x, merge_df$FOT.y)
cor.g <- ggpairs(
  log10(merge_df[,c(4,7)]),
  columnLabels = c('QE-Plus-Global', 'QE-HF-Global'),
  lower=list(continuous=my_fn)
)


# High
merge_df.high = merge_df[which(
  merge_df$ExprLab.x=='High' & merge_df$ExprLab.y=='High'
),]
cor(merge_df.high$FOT.x, merge_df.high$FOT.y)
cor.h <- ggpairs(
  log10(merge_df.high[,c(4,7)]),
  columnLabels = c('QE-Plus-High', 'QE-HF-High')
)

# Medium
merge_df.medium = merge_df[which(
  merge_df$ExprLab.x=='Medium' & merge_df$ExprLab.y=='Medium'
),]
cor(merge_df.medium$FOT.x, merge_df.medium$FOT.y)
cor.m <- ggpairs(
  log10( merge_df.medium[,c(4,7)]),
  columnLabels = c('QE-Plus-Medium', 'QE-HF-Medium')
)

# Low
merge_df.low = merge_df[which(
  merge_df$ExprLab.x=='Low' & merge_df$ExprLab.y=='Low'
),]
cor(merge_df.low$FOT.x, merge_df.low$FOT.y)
cor.l <- ggpairs(
  log10( merge_df.low[,c(4,7)]),
  columnLabels = c('QE-Plus-Low', 'QE-HF-Low')
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
  p <- p + xlab("log10(FOT)/QE-Plus") + ylab("log10(FOT)/QE-HF") 
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
h_index = which(merge_df$ExprLab.x=='High' & merge_df$ExprLab.y=='High')
m_index = which(merge_df$ExprLab.x=='Medium' & merge_df$ExprLab.y=='Medium')
l_index = which(merge_df$ExprLab.x=='Low' & merge_df$ExprLab.y=='Low')
qeplus.px.fot.common = qeplus.px.fot[match(common_genes,qeplus.px.fot$GeneSymbol),]
qehf.px.fot.common = qehf.px.fot[match(common_genes,qehf.px.fot$GeneSymbol),]

get_cv_list <- function(df){
  # df = qeplus.px.fot.common
  df.row_cv = apply(df[,-1], 1, function(x){
    sd(x)/mean(x)
  })
  df.row_cv = round(df.row_cv, 2)
  return(df.row_cv)
}

QE_Plus.g.cv = get_cv_list(qeplus.px.fot.common)
QE_Plus.h.cv = get_cv_list(qeplus.px.fot.common[h_index,])
QE_Plus.m.cv = get_cv_list(qeplus.px.fot.common[m_index,])
QE_Plus.l.cv = get_cv_list(qeplus.px.fot.common[l_index,])

QE_HF.g.cv = get_cv_list(qehf.px.fot.common)
QE_HF.h.cv = get_cv_list(qehf.px.fot.common[h_index,])
QE_HF.m.cv = get_cv_list(qehf.px.fot.common[m_index,])
QE_HF.l.cv = get_cv_list(qehf.px.fot.common[l_index,])

QE_HFX_QE_HF.g.cv = get_cv_list(
  data.frame(qeplus.px.fot.common, qehf.px.fot.common[,-1])
)
QE_HFX_QE_HF.h.cv = get_cv_list(
  data.frame(qeplus.px.fot.common, qehf.px.fot.common[,-1])[h_index,]
)
QE_HFX_QE_HF.m.cv = get_cv_list(
  data.frame(qeplus.px.fot.common, qehf.px.fot.common[,-1])[m_index,]
)
QE_HFX_QE_HF.l.cv = get_cv_list(
  data.frame(qeplus.px.fot.common, qehf.px.fot.common[,-1])[l_index,]
)

bp <- boxplot(
  QE_Plus.l.cv, QE_Plus.m.cv, QE_Plus.h.cv, QE_Plus.g.cv,
  QE_HF.l.cv, QE_HF.m.cv, QE_HF.h.cv, QE_HF.g.cv,
  QE_HFX_QE_HF.l.cv, QE_HFX_QE_HF.m.cv, QE_HFX_QE_HF.h.cv, QE_HFX_QE_HF.g.cv
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
    QE_Plus.l.cv, QE_Plus.m.cv, QE_Plus.h.cv, QE_Plus.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25), boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  boxplot(
    QE_HF.l.cv, QE_HF.m.cv, QE_HF.h.cv, QE_HF.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25)+4.5, boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  boxplot(
    QE_HFX_QE_HF.l.cv, QE_HFX_QE_HF.m.cv, QE_HFX_QE_HF.h.cv, QE_HFX_QE_HF.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25)+4.5*2, boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  abline(v = 4, lty = 3, lwd = 2, col = 'grey')
  abline(v = 8.5, lty = 3, lwd = 2, col = 'grey')
  axis(
    1, at = c(1.75, 6.25, 10.75), labels = c('QE-Plus', 'QE-HF', 'QE-Plus&QE-HF')
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
# meta_data.NPB.QEPlus 
# meta_data.NPB.QEHF

# D5
QE_Plus.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEPlus$Exp_code)
QE_HF.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEHF$Exp_code)

get_deps <- function(df, group, min_fc=2){
  # df = QE_Plus.fot
  # group = meta_data.NPB.QEPlus$UPerson
  
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

intersect_genes = intersect(as.vector(QE_Plus.fot$GeneSymbol), as.vector(QE_HF.fot$GeneSymbol))
min_fc = 2
QE_Plus.fot.deps_list = get_deps(QE_Plus.fot[match(intersect_genes, QE_Plus.fot$GeneSymbol),], meta_data.NPB.QEPlus$UPerson, min_fc=2)
QE_HF.fot.deps_list = get_deps(QE_HF.fot[match(intersect_genes, QE_HF.fot$GeneSymbol),], meta_data.NPB.QEHF$UPerson, min_fc=2)

get_consensus_deps = function(df1, df2, min_fc = 2){
  df1 = QE_Plus.fot.deps_list[[1]]
  df2 = QE_HF.fot.deps_list[[1]]
  
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
    s_lab = paste('QE-Plus', '-', s_lab, sep = '')
    t_lab = paste('QE-HF', '-', t_lab, sep = '')
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
    links_df$source, levels = c('QE-Plus-DOWN', 'QE-Plus-NA', 'QE-Plus-UP')
  )
  links_df$target = factor(
    links_df$target, levels = c('QE-HF-UP', 'QE-HF-DOWN', 'QE-HF-NA')
  )
  
  links = links_df[order(links_df$source),]
  # links = links[c(3,1,2, 6,4,5, 8,7),]
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )
  # nodes$name = factor(nodes$name, levels = c('QE-Plus-UP', 'QE-Plus-DOWN', 'QE-Plus-NA', 'QE-HF-UP', 'QE-HF-DOWN', 'QE-HF-NA'))
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Make the Network
  # '#33C860', # 绿色
  # '#81B0FF', # 蓝色
  # '#F9918A'  # 红色
  my_color <- 'd3.scaleOrdinal() .domain(["QE-Plus-UP", "QE-Plus-DOWN", "QE-Plus-NA", "QE-HF-UP", "QE-HF-DOWN", "QE-HF-NA"]) .range(["#F9918A",  "#33C860", "grey", "#F9918A",  "#33C860", "grey"])'
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



