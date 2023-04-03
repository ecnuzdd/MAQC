# fot_data
# meta_data

meta_data.NPB.Fusion = meta_data[which(
  meta_data$Labcode == 'NPB' 
  & meta_data$SInstrument == 'Fusion'
  & meta_data$Purpose == 'Replications'
),]



meta_data.NPB.Lumos = meta_data[which(
  meta_data$Labcode == 'NPB' 
  & meta_data$SInstrument == 'Lumos'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]

# D5
px_list = c('D5', 'D6', 'F7', 'M8')
px = px_list[1]
fusion.px.index = which(meta_data.NPB.Fusion$UPerson == px)
fusion.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Fusion$Exp_code[fusion.px.index])
delele.index = which(fusion.px.fot$GeneSymbol ==	'ACTA1') # contamination
if(length(which(delele.index>0))){
  print(delele.index)
  fusion.px.fot = fusion.px.fot[-delele.index,] # 峰度排名在单个实验排名第一，但只被鉴定到一次
}

lumos.px.index = which(meta_data.NPB.Lumos$UPerson == px)
lumos.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[lumos.px.index])

#1 identification
fusion.px.col_gps = apply(fusion.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(fusion.px.fot[,1])))
lumos.px.col_gps = apply(lumos.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(lumos.px.fot[,1])))


fusion.px.id = as.vector(unlist(lapply(fusion.px.col_gps, length)))
lumos.px.id = as.vector(unlist(lapply(lumos.px.col_gps, length)))
barplot_df = data.frame(
  ID =  c(fusion.px.id, nrow(fusion.px.fot), lumos.px.id, nrow(lumos.px.fot)),
  Exptype = rep(c('R1', 'R2', 'R3', 'Union'), 2),
  Instrument = rep(c('Fusion', 'Lumos'), each=4)
)
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
  ylim = c(0, 5000),
  font.label = list(color = "black", size = 10, vjust = -1)
)

# venn
library(VennDiagram)
library(tidyverse)
alpha_per = 0.3
jco_col = pal_jco("default", alpha = 0.5)(10)[c(1)]
venn.diagram(
  x = list(
    R1 = fusion.px.col_gps[[1]],
    R2 = fusion.px.col_gps[[2]],
    R3 = fusion.px.col_gps[[3]]
  ),
  filename = 'temp_figures/03/NPB_Fusion_Lumos/03_venn_Fusion.png',
  main = 'Fusion',
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
    R1 = lumos.px.col_gps[[1]],
    R2 = lumos.px.col_gps[[2]],
    R3 = lumos.px.col_gps[[3]]
  ),
  filename = 'temp_figures/03/NPB_Fusion_Lumos/03_venn_Lumos.png',
  main = 'Lumos',
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

#Fusion-Lumos
jco_col = pal_jco("default", alpha = 0.5)(10)[c(1,2)]
venn.diagram(
  x = list(
    Lumos = as.vector(lumos.px.fot$GeneSymbol),
    Fusion = as.vector(fusion.px.fot$GeneSymbol)
  ),
  filename = 'temp_figures/03/NPB_Fusion_Lumos/03_venn_Fusion_Lumos_per.png',
  main = 'Fusion vs Lumos',
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
fusion.row_id = apply(fusion.px.fot[,-1], 1, function(x){length(which(x>0))})
lumos.row_id = apply(lumos.px.fot[,-1], 1, function(x){length(which(x>0))})













#3 intensity distribution
get_expr_label <- function(x){
  # x = fusion.row_mean
  x.quantile = as.vector(quantile(x, probs = c(1/3, 2/3)))
  x.l_index = which(x<x.quantile[1])
  x.h_index = which(x>x.quantile[2])
  x.expr_label = rep('Medium', length(x))
  x.expr_label[x.l_index] = 'Low'
  x.expr_label[x.h_index] = 'High'
  x.expr_label = factor(x.expr_label, levels = c('Low', 'Medium', 'High'))
  return(x.expr_label)
}

fusion.row_mean = rowMeans(fusion.px.fot[,-1])
fusion.row_expr_label = get_expr_label(fusion.row_mean)

lumos.row_mean = rowMeans(lumos.px.fot[,-1])
lumos.row_expr_label = get_expr_label(lumos.row_mean)

fusion.df = data.frame(
  GeneSymbol = fusion.px.fot$GeneSymbol,
  RowID = fusion.row_id,
  ExprLab = fusion.row_expr_label
)

lumos.df = data.frame(
  GeneSymbol = lumos.px.fot$GeneSymbol,
  RowID = lumos.row_id,
  ExprLab = lumos.row_expr_label
)

table(fusion.df$RowID, fusion.df$ExprLab)
table(lumos.df$RowID, lumos.df$ExprLab)
library(ggforce)
df = fusion.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
PieDonut_zdd(
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(1,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = 'Fusion (Total=4202)'
)


df = lumos.df
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
  title = 'Lumos (Total=4425)'
)






#2 correlation
#fusion
if(T){
  df.ref = fusion.df
  df = fusion.px.fot
  
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
  pdf('03_pairwise_cor_Fusion.pdf', height = 5, width = 5)
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
  #lumos
  df.ref = lumos.df
  df = lumos.px.fot
  
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
  
  pdf('03_pairwise_cor_Lumos.pdf', height = 5, width = 5)
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


# fusion.px.col_gps, fusion.df
fusion.repro = get_px_reproducibility(fusion.px.col_gps, ref_df = NULL, ref_lab = NULL)
fusion.repro.h = get_px_reproducibility(fusion.px.col_gps, ref_df = fusion.df, ref_lab = 'High')
fusion.repro.m = get_px_reproducibility(fusion.px.col_gps, ref_df = fusion.df, ref_lab = 'Medium')
fusion.repro.l = get_px_reproducibility(fusion.px.col_gps, ref_df = fusion.df, ref_lab = 'Low')
df.fusion = data.frame(
  Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
  Reproducibity = c(fusion.repro.l, fusion.repro.m, fusion.repro.h, fusion.repro)
)
df.fusion$Intensity = factor(df.fusion$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))


# lumos.px.col_gps, lumos.df
lumos.repro = get_px_reproducibility(lumos.px.col_gps, ref_df = NULL, ref_lab = NULL)
lumos.repro.h = get_px_reproducibility(lumos.px.col_gps, ref_df = lumos.df, ref_lab = 'High')
lumos.repro.m = get_px_reproducibility(lumos.px.col_gps, ref_df = lumos.df, ref_lab = 'Medium')
lumos.repro.l = get_px_reproducibility(lumos.px.col_gps, ref_df = lumos.df, ref_lab = 'Low')
df.lumos = data.frame(
  Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
  Reproducibity = c(lumos.repro.l, lumos.repro.m, lumos.repro.h, lumos.repro)
)
df.lumos$Intensity = factor(df.lumos$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))

df = rbind.data.frame(
  df.fusion, 
  df.lumos
)
df$Instrument = rep(
  c('Fusion', 'Lumos'), each = 24
)
df$Instrument = factor(df$Instrument, levels = c('Fusion', 'Lumos'))

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
# fusion.px.fot
# lumos.px.fot
if(T){
  fusion.px.fot.row_id = apply(fusion.px.fot[,-1], 1, function(x){length(which(x>0))})
  lumos.px.fot.row_id = apply(lumos.px.fot[,-1], 1, function(x){length(which(x>0))})
  px.fot = merge(
    fusion.px.fot[fusion.px.fot.row_id>=1,], 
    lumos.px.fot[lumos.px.fot.row_id>=1,], 
    by = 'GeneSymbol', all = T
  )
  px.fot.mat = as.matrix(px.fot[,-1])
  colnames(px.fot.mat) = c(
    paste('Fusion', c('R1', 'R2', 'R3'), sep = '_'),
    paste('Lumos', c('R1', 'R2', 'R3'), sep = '_')
  )
  na.index = which(is.na(px.fot.mat))
  px.fot.mat[na.index] = 0
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
    limits = c(0.8, 1), breaks = seq(0.8,1, 0.05),
    guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)
  ) + labs(fill = "PCorr")
  p1
  
  # col<- colorRampPalette(c("blue", "white", "red"))(20)
  # heatmap(x = cormat, col = col, symm = TRUE)
}














# High
fusion.h.gps = fusion.df$GeneSymbol[fusion.df$ExprLab=='High']
lumos.h.gps = lumos.df$GeneSymbol[lumos.df$ExprLab=='High']
length(intersect(fusion.h.gps, lumos.h.gps))/length(union(fusion.h.gps, lumos.h.gps))

# Medium
fusion.m.gps = fusion.df$GeneSymbol[fusion.df$ExprLab=='Medium']
lumos.m.gps = lumos.df$GeneSymbol[lumos.df$ExprLab=='Medium']
length(intersect(fusion.m.gps, lumos.m.gps))/length(union(fusion.m.gps, lumos.m.gps))

# Low
fusion.l.gps = fusion.df$GeneSymbol[fusion.df$ExprLab=='Low']
lumos.l.gps = lumos.df$GeneSymbol[lumos.df$ExprLab=='Low']
length(intersect(fusion.l.gps, lumos.l.gps))/length(union(fusion.l.gps, lumos.l.gps))



# Sankey
sankey_df = merge(fusion.df, lumos.df, by = 'GeneSymbol', all = T)
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
  s_lab = paste('Fusion', '-', s_lab, sep = '')
  t_lab = paste('Lumos', '-', t_lab, sep = '')
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
  links_df$source, levels = c('Fusion-High', 'Fusion-Medium', 'Fusion-Low', 'Fusion-NA')
)
links_df$target = factor(
  links_df$target, levels = c('Lumos-High', 'Lumos-Medium', 'Lumos-Low', 'Lumos-NA')
)

links = links_df[order(links_df$source),]
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
my_color <- 'd3.scaleOrdinal() .domain(["Fusion-High", "Fusion-Medium", "Fusion-Low", "Fusion-NA", "Lumos-High", "Lumos-Medium", "Lumos-Low", "Lumos-NA"]) .range(["#F9918A", "#81B0FF", "#33C860", "grey", "#F9918A", "#81B0FF", "#33C860", "grey"])'
p <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", 
  NodeID = "name", 
  fontSize = 15,
  colourScale = my_color,
  LinkGroup="group",
  sinksRight = FALSE
)
p

# 相关性
fusion.df$FOT = apply(fusion.px.fot[,-1], 1, mean)
lumos.df$FOT = apply(lumos.px.fot[,-1], 1, mean)
merge_df = merge(fusion.df, lumos.df, by = 'GeneSymbol', all = T)
merge_df = na.omit(merge_df)

# Global
cor(merge_df$FOT.x, merge_df$FOT.y)
cor.g <- ggpairs(
  log10(merge_df[,c(4,7)]),
  columnLabels = c('Fusion-Global', 'Lumos-Global'),
  lower=list(continuous=my_fn)
)


# High
merge_df.high = merge_df[which(
  merge_df$ExprLab.x=='High' & merge_df$ExprLab.y=='High'
),]
cor(merge_df.high$FOT.x, merge_df.high$FOT.y)
cor.h <- ggpairs(
  log10(merge_df.high[,c(4,7)]),
  columnLabels = c('Fusion-High', 'Lumos-High')
)

# Medium
merge_df.medium = merge_df[which(
  merge_df$ExprLab.x=='Medium' & merge_df$ExprLab.y=='Medium'
),]
cor(merge_df.medium$FOT.x, merge_df.medium$FOT.y)
cor.m <- ggpairs(
  log10( merge_df.medium[,c(4,7)]),
  columnLabels = c('Fusion-Medium', 'Lumos-Medium')
)

# Low
merge_df.low = merge_df[which(
  merge_df$ExprLab.x=='Low' & merge_df$ExprLab.y=='Low'
),]
cor(merge_df.low$FOT.x, merge_df.low$FOT.y)
cor.l <- ggpairs(
  log10( merge_df.low[,c(4,7)]),
  columnLabels = c('Fusion-Low', 'Lumos-Low')
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
  p <- p + xlab("log10(FOT)/Fusion") + ylab("log10(FOT)/Lumos") 
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

ggarrange(
  p.l, p.m, p.h, p.g, 
  nrow = 2, ncol = 2,
  labels = c('A', 'B', 'C', 'D')
)


# CV
common_genes = as.vector(merge_df$GeneSymbol)
h_index = which(merge_df$ExprLab.x=='High' & merge_df$ExprLab.y=='High')
m_index = which(merge_df$ExprLab.x=='Medium' & merge_df$ExprLab.y=='Medium')
l_index = which(merge_df$ExprLab.x=='Low' & merge_df$ExprLab.y=='Low')
fusion.px.fot.common = fusion.px.fot[match(common_genes,fusion.px.fot$GeneSymbol),]
lumos.px.fot.common = lumos.px.fot[match(common_genes,lumos.px.fot$GeneSymbol),]

get_cv_list <- function(df){
  # df = fusion.px.fot.common
  df.row_cv = apply(df[,-1], 1, function(x){
    sd(x)/mean(x)
  })
  df.row_cv = round(df.row_cv, 2)
  return(df.row_cv)
}

fusion.g.cv = get_cv_list(fusion.px.fot.common)
fusion.h.cv = get_cv_list(fusion.px.fot.common[h_index,])
fusion.m.cv = get_cv_list(fusion.px.fot.common[m_index,])
fusion.l.cv = get_cv_list(fusion.px.fot.common[l_index,])

lumos.g.cv = get_cv_list(lumos.px.fot.common)
lumos.h.cv = get_cv_list(lumos.px.fot.common[h_index,])
lumos.m.cv = get_cv_list(lumos.px.fot.common[m_index,])
lumos.l.cv = get_cv_list(lumos.px.fot.common[l_index,])

fusion_lumos.g.cv = get_cv_list(
  data.frame(fusion.px.fot.common, lumos.px.fot.common[,-1])
)
fusion_lumos.h.cv = get_cv_list(
  data.frame(fusion.px.fot.common, lumos.px.fot.common[,-1])[h_index,]
)
fusion_lumos.m.cv = get_cv_list(
  data.frame(fusion.px.fot.common, lumos.px.fot.common[,-1])[m_index,]
)
fusion_lumos.l.cv = get_cv_list(
  data.frame(fusion.px.fot.common, lumos.px.fot.common[,-1])[l_index,]
)

bp <- boxplot(
  fusion.l.cv, fusion.m.cv, fusion.h.cv, fusion.g.cv,
  lumos.l.cv, lumos.m.cv, lumos.h.cv, lumos.g.cv,
  fusion_lumos.l.cv, fusion_lumos.m.cv, fusion_lumos.h.cv, fusion_lumos.g.cv
)

# '#33C860', # 绿色
# '#81B0FF', # 蓝色
# '#F9918A'  # 红色

if(T){
  bp <- boxplot(
    NA, 
    xlim = c(0, 12.5), ylim = c(0, 3), ylab = 'Coefficient of Variation',
    xaxt = 'n', yaxt = 'n'
  )
  boxplot(
    fusion.l.cv, fusion.m.cv, fusion.h.cv, fusion.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25), boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  boxplot(
    lumos.l.cv, lumos.m.cv, lumos.h.cv, lumos.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25)+4.5, boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  boxplot(
    fusion_lumos.l.cv, fusion_lumos.m.cv, fusion_lumos.h.cv, fusion_lumos.g.cv,
    at = c(0.25, 1.25, 2.25, 3.25)+4.5*2, boxwex = 0.5,
    add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
  )
  
  abline(v = 4, lty = 3, lwd = 2, col = 'grey')
  abline(v = 8.5, lty = 3, lwd = 2, col = 'grey')
  axis(
    1, at = c(1.75, 6.25, 10.75), labels = c('Fusion', 'Lumos', 'Fusion&Lumos')
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
# meta_data.NPB.Fusion 
# meta_data.NPB.Lumos

# D5
fusion.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Fusion$Exp_code)
lumos.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code)

get_deps <- function(df, group, min_fc=2){
  # df = fusion.fot
  # group = meta_data.NPB.Fusion$UPerson
  
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

fusion.fot.deps_list = get_deps(fusion.fot, meta_data.NPB.Fusion$UPerson, min_fc=2)
lumos.fot.deps_list = get_deps(lumos.fot, meta_data.NPB.Lumos$UPerson, min_fc=2)

get_consensus_deps = function(df1, df2, min_fc = 2){
  df1 = fusion.fot.deps_list[[1]]
  df2 = lumos.fot.deps_list[[1]]
  c_gps.h = intersect(
    df1$GeneSymbol[df1$FC>min_fc], df2$GeneSymbol[df2$FC>1/min_fc]
  )
  c_gps.l = intersect(
    df1$GeneSymbol[df1$FC<min_fc], df2$GeneSymbol[df2$FC<1/min_fc]
  )
  c_gps = c(c_gps.h, c_gps.l)
  
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
  sankey_df = sankey_df[,c(1,3,5)]
  
  a = sankey_df[which(
    (sankey_df$Ctrl_id.x>1 | sankey_df$Case_id.x>1)
    | (sankey_df$Ctrl_id.y>1 | sankey_df$Case_id.y>1)
  ),]
  
  a.index = which(
    a$Label.x == a$Label.y
  )

  print(length(a.index)/nrow(a))
  
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
    s_lab = paste('Fusion', '-', s_lab, sep = '')
    t_lab = paste('Lumos', '-', t_lab, sep = '')
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
  
  # links_df$source = factor(
  #   links_df$source, levels = c('Fusion-UP', 'Fusion-DOWN', 'Fusion-NA')
  # )
  # links_df$target = factor(
  #   links_df$target, levels = c('Lumos-UP', 'Lumos-DOWN', 'Lumos-NA')
  # )
  
  links = links_df[order(links_df$source),]
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )
  # nodes$name = factor(nodes$name, levels = c('Fusion-UP', 'Fusion-DOWN', 'Fusion-NA', 'Lumos-UP', 'Lumos-DOWN', 'Lumos-NA'))
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Make the Network
  # '#33C860', # 绿色
  # '#81B0FF', # 蓝色
  # '#F9918A'  # 红色
  my_color <- 'd3.scaleOrdinal() .domain(["Fusion-UP", "Fusion-DOWN", "Fusion-NA", "Lumos-UP", "Lumos-DOWN", "Lumos-NA"]) .range(["#F9918A",  "#33C860", "grey", "#F9918A",  "#33C860", "grey"])'
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







# PCA
fusion.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Fusion$Exp_code)
lumos.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code)

# fusion.fot
fusion.fot.row_mean = rowMeans(fusion.fot[,-1])
fusion.fot.row_expr_label = get_expr_label(fusion.fot.row_mean)
fusion.fot.df = data.frame(
  GeneSymbol = fusion.fot$GeneSymbol,
  RowID = apply(fusion.fot[,-1], 1, function(x){length(which(x>0))}),
  MeanFOT = fusion.fot.row_mean,
  ExprLab = fusion.fot.row_expr_label
)

# lumos.fot
lumos.fot.row_mean = rowMeans(lumos.fot[,-1])
lumos.fot.row_expr_label = get_expr_label(lumos.fot.row_mean)
lumos.fot.df = data.frame(
  GeneSymbol = lumos.fot$GeneSymbol,
  RowID = apply(lumos.fot[,-1], 1, function(x){length(which(x>0))}),
  MeanFOT = lumos.fot.row_mean,
  ExprLab = lumos.fot.row_expr_label
)

# Single PCA
# Fusion
if(T){
  fot_df = fusion.fot
  fot_df.mat = as.matrix(fot_df[,-1])
  fot_df.mat.scale = t(scale(t(fot_df.mat)))
  fot_df.ref = fusion.fot.df
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
    main = paste('Fusion; Low-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
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
    main = paste('Fusion; Medium-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
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
    main = paste('Fusion; High-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
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
    main = paste('Fusion; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
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


# Lumos
if(T){
  fot_df = lumos.fot
  fot_df.ref = lumos.fot.df
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
    main = paste('Lumos; Low-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
    sub = 'Lumos-Lumos'
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
    main = paste('Lumos; Medium-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
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
    main = paste('Lumos; High-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
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
    main = paste('Lumos; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
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




fusion.fot_df = fusion.fot
fusion.fot_df.ref = fusion.fot.df
fusion.fot_df.mat = as.matrix(fusion.fot_df[,-1])
fusion.fot_df.mat.scale = t(scale(t(fusion.fot_df.mat)))
fusion.fot_group = rep(p_flags, each = 3)

lumos.fot_df = lumos.fot
lumos.fot_df.ref = lumos.fot.df
lumos.fot_df.mat = as.matrix(lumos.fot_df[,-1])
lumos.fot_df.mat.scale = t(scale(t(lumos.fot_df.mat)))
lumos.fot_group = rep(p_flags, each = 3)

get_merge_df <- function(df1, df2, replace_na = T){
  # df1 = fusion.fot_df
  # df2 = lumos.fot_df
  merge_df = merge(df1, df2, by = 'GeneSymbol', all = T)
  merge_df.mat = as.matrix(merge_df[,-1])
  if(replace_na){
    na.index = which(is.na(merge_df.mat))
    merge_df.mat[na.index] = 0
  }
  merge_df = data.frame(
    GeneSymbol = merge_df$GeneSymbol,
    merge_df.mat
  )
  return(merge_df)
}


if(T){ # 8*8
  # Global
  merge_df = get_merge_df(fusion.fot_df, lumos.fot_df)
  merge_df.mat = as.matrix(merge_df[,-1])
  merge_df.mat.scale = t(scale(t(merge_df.mat)))
  merge_df_group = rep(
    rep(p_flags, each = 3), 2
  )
  merge_df_type = rep(c(16,17), each = 12)
  
  mat = merge_df.mat.scale #[fot_df.ref$ExprLab=='Low',]
  snr = round(snr_function(mat, merge_df_group),2)
  colors = rep(
    rep(p_colors, each = 3), 2
  )
  pca <- prcomp(((t(mat))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste('Fusion_Lumos; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
  plot(
    pca_predict.2d, 
    t = 'n', 
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d, pch = merge_df_type, col = colors)
}




if(T){ # 8*8
  # Global
  fusion.fot_df.scale = data.frame(
    GeneSymbol = fusion.fot_df$GeneSymbol,
    t(scale(t(fusion.fot_df[,-1])))
  )
  lumos.fot_df.scale = data.frame(
    GeneSymbol = lumos.fot_df$GeneSymbol,
    t(scale(t(lumos.fot_df[,-1])))
  )
  merge_df = get_merge_df(fusion.fot_df.scale, lumos.fot_df.scale, replace_na = F)
  merge_df = na.omit(merge_df)
  merge_df.mat = as.matrix(merge_df[,-1])
  merge_df.mat.scale = merge_df.mat
  merge_df_group = rep(
    rep(p_flags, each = 3), 2
  )
  merge_df_type = rep(c(16,17), each = 12)
  
  mat = merge_df.mat.scale #[fot_df.ref$ExprLab=='Low',]
  snr = round(snr_function(mat, merge_df_group),2)
  colors = rep(
    rep(p_colors, each = 3), 2
  )
  pca <- prcomp(((t(mat))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste('Fusion_Lumos; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
  plot(
    pca_predict.2d, 
    t = 'n', 
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d, pch = merge_df_type, col = colors)
}


#### CV ####
if(T){
  
  
  meta_data.NPB.Fusion = meta_data[which(
    meta_data$Labcode == 'NPB' 
    & meta_data$SInstrument == 'Fusion'
    & meta_data$Purpose == 'Replications'
  ),]
  
  
  
  meta_data.NPB.Lumos = meta_data[which(
    meta_data$Labcode == 'NPB' 
    & meta_data$SInstrument == 'Lumos'
    & meta_data$Purpose == 'Replications'
    & meta_data$Batch == 1
  ),]
  
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px = px_list[1]
  fusion.px.index = which(meta_data.NPB.Fusion$UPerson == px)
  fusion.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Fusion$Exp_code[fusion.px.index])
  # delele.index = which(fusion.px.fot$GeneSymbol ==	'ACTA1') # contamination
  # if(length(which(delele.index>0))){
  #   print(delele.index)
  #   fusion.px.fot = fusion.px.fot[-delele.index,] # 峰度排名在单个实验排名第一，但只被鉴定到一次
  # }
  # 
  lumos.px.index = which(meta_data.NPB.Lumos$UPerson == px)
  lumos.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[lumos.px.index])
  
  
  fusion.px.fot.cv_list = get_cv_list_set(fusion.px.fot, min_id = 3)
  lumos.px.fot.cv_list = get_cv_list_set(lumos.px.fot, min_id = 3)
  
  
  
  
  # '#33C860', # 绿色
  # '#81B0FF', # 蓝色
  # '#F9918A'  # 红色
  
  if(T){
    bp <- boxplot(
      NA, 
      xlim = c(0, 8), ylim = c(0, 2), ylab = 'Coefficient of Variation',
      xaxt = 'n', yaxt = 'n'
    )
    boxplot(
      fusion.px.fot.cv_list$l_cv, fusion.px.fot.cv_list$m_cv,
      fusion.px.fot.cv_list$h_cv, fusion.px.fot.cv_list$g_cv,
      at = c(0.25, 1.25, 2.25, 3.25), boxwex = 0.5,
      add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
    )
    
    boxplot(
      lumos.px.fot.cv_list$l_cv, lumos.px.fot.cv_list$m_cv,
      lumos.px.fot.cv_list$h_cv, lumos.px.fot.cv_list$g_cv,
      at = c(0.25, 1.25, 2.25, 3.25)+4.5, boxwex = 0.5,
      add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
    )
    
    
    
    
    abline(v = c(4, 8.5, 13, 17.5)[1], lty = 3, lwd = 2, col = 'grey')
    axis(
      1, at = c(1.75, 6.25, 10.75, 15.25, 19.75)[1:2], 
      labels = c('Fusion', 'Lumos', 'G3', 'G4', 'G5')[1:2]
    )
    
    axis(
      2, at = seq(0, 3, 0.5), 
      labels = c('0', '50%', '100%', '150%', '200%', '250%', '300%'), 
      las = 1, cex.axis = 0.8
    )
    axis(
      2, at = c(0.15, 0.2, 0.3), 
      labels = c('15%', '20%', '30%'), 
      tcl = -0.3, cex.axis = 0.6, las = 1, col.axis = 'red', col = 'red'
    )
    abline(h = c(0.15, 0.2, 0.3), lty = 3, lwd = 2, col = 'red')
  }
  
}








if(F){
  get_probs_table <- function(df){
    # qe_qeplus.df
    # qeplus_qehf.df
    # qehf_qehfx.df
    # qehfx_e480.df
    df = df[,c('source', 'target', 'value')]
    del.index = apply(df, 1, function(x){
      x = as.vector(unlist(x))
      if(grepl('NA', x[1])){
        return(T)
      }else{
        return(F)
      }
    })
    df1 = df[-which(del.index),]
    df1.mat = matrix(as.numeric(as.vector(df1$value)), nrow = 3, ncol = 4, byrow = T)
    rownames(df1.mat) = unique(df1$source)
    colnames(df1.mat) = unique(df1$target)
    df1.mat.rowsum = apply(df1.mat, 1, sum)
    df1.mat.probs_table = df1.mat/df1.mat.rowsum
    df1.mat.probs_table = round(df1.mat.probs_table, 2)
    return(df1.mat.probs_table)
  }
  
  fusion_lumos.df = links
  fusion_lumos.df.table = get_probs_table(fusion_lumos.df)
  fusion_lumos.df.table = fusion_lumos.df.table[,c(1,3,2,4)]
  output_path = normalizePath(path = file.path(
    BASE_DIR, 'sankey','Fusion-Lumos-sankey-probs.xlsx'
  ),mustWork = F)
  
  write.xlsx(
    fusion_lumos.df.table, output_path, sheetName = 'fusion_lumos'
  )
 
  
  
  
  
  # PCA
  fusion.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Fusion$Exp_code)
  lumos.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code)
  
  fusion.fot.ref = get_fot_ref_df(fusion.fot)
  lumos.fot.ref = get_fot_ref_df(lumos.fot)

  pdf(
    normalizePath(
      file.path(BASE_DIR, 'sankey', '20200825', 'PCA_fusion_lumos.pdf'), mustWork = F
    ),
    height = 8, width = 8
  )
  g1.snr_list = plot_single_pca(fusion.fot, fusion.fot.ref, main_lab = 'Fusion', p_flags, p_colors)
  g2.snr_list = plot_single_pca(lumos.fot, lumos.fot.ref, main_lab = 'Lumos', p_flags, p_colors)

  dev.off()
  
}
