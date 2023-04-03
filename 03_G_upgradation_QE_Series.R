# fot_data
# meta_data

meta_data.NPS.QE = meta_data[which(
  meta_data$Labcode == 'NPS' 
  & meta_data$SInstrument == 'QE'
  & meta_data$Purpose == 'Replications'
),]

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


meta_data.TMO.QEHFX = meta_data[which(
  meta_data$Labcode == 'TMO' 
  & meta_data$SInstrument == 'QE-HFX'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]


meta_data.TMO.E480 = meta_data[which(
  meta_data$Labcode == 'TMO' 
  & meta_data$SInstrument == 'Exploris480'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]








# D5
px_list = c('D5', 'D6', 'F7', 'M8')
px = px_list[1]

# QE
qe.px.index = which(meta_data.NPS.QE$UPerson == px)
qe.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPS.QE$Exp_code[qe.px.index])

# QE-Plus
qeplus.px.index = which(meta_data.NPB.QEPlus$UPerson == px)
qeplus.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEPlus$Exp_code[qeplus.px.index])

# QE-HF
qehf.px.index = which(meta_data.NPB.QEHF$UPerson == px)
qehf.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEHF$Exp_code[qehf.px.index])

# QE-HFX
qehfx.px.index = which(meta_data.TMO.QEHFX$UPerson == px)
qehfx.px.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.QEHFX$Exp_code[qehfx.px.index])

# Exploris480
e480.px.index = which(meta_data.TMO.E480$UPerson == px)
e480.px.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.E480$Exp_code[e480.px.index])


#1 identification
qe.px.col_gps = apply(qe.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(qe.px.fot[,1])))
qeplus.px.col_gps = apply(qeplus.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(qeplus.px.fot[,1])))
qehf.px.col_gps = apply(qehf.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(qehf.px.fot[,1])))
qehfx.px.col_gps = apply(qehfx.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(qehfx.px.fot[,1])))
e480.px.col_gps = apply(e480.px.fot[,-1], 2, function(x, y){
  y[which(x>0)]
}, y = as.vector(unlist(e480.px.fot[,1])))

qe.px.id = as.vector(unlist(lapply(qe.px.col_gps, length)))
qeplus.px.id = as.vector(unlist(lapply(qeplus.px.col_gps, length)))
qehf.px.id = as.vector(unlist(lapply(qehf.px.col_gps, length)))
qehfx.px.id = as.vector(unlist(lapply(qehfx.px.col_gps, length)))
e480.px.id = as.vector(unlist(lapply(e480.px.col_gps, length)))

barplot_df = data.frame(
  ID =  c(
    qe.px.id, nrow(qe.px.fot),
    qeplus.px.id, nrow(qeplus.px.fot), qehf.px.id, nrow(qehf.px.fot),
    qehfx.px.id, nrow(qehfx.px.fot), e480.px.id, nrow(e480.px.fot)
  ),
  Exptype = rep(c('R1', 'R2', 'R3', 'Union'), 5),
  Instrument = rep(c('QE', 'QE-Plus', 'QE-HF', 'QE-HFX', 'Exploris480'), each=4)
)
barplot_df$Class = paste(barplot_df$Instrument, barplot_df$Exptype, sep = '_')
barplot_df$Class = factor(barplot_df$Class, levels = barplot_df$Class)
barplot_df$Instrument = factor(barplot_df$Instrument, levels = c('QE', 'QE-Plus', 'QE-HF', 'QE-HFX', 'Exploris480'))

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

jco_col = pal_jco("default", alpha = 0.5)(10)[1:5]
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
  add.params = list(color = rep(c(jco_col), each = 4), size = 2), # Change segment color and size,
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
pdf(file="temp_figures/03/QE_series/03_venn_QE_S.pdf")

temp <- venn.diagram(
  x = list(
    R1 = qe.px.col_gps[[1]],
    R2 = qe.px.col_gps[[2]],
    R3 = qe.px.col_gps[[3]]
  ),
  # filename = 'temp_figures/03/QE_series/03_venn_QE_S.pdf',
  filename = NULL,
  main = 'qe',
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
grid.draw(temp)
dev.off()
jco_col = pal_jco("default", alpha = 0.5)(10)[c(2)]
venn.diagram(
  x = list(
    R1 = qeplus.px.col_gps[[1]],
    R2 = qeplus.px.col_gps[[2]],
    R3 = qeplus.px.col_gps[[3]]
  ),
  filename = 'Venn/03_venn_qeplus.png',
  main = 'qeplus',
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

#qe-qeplus
library(grDevices)
pdf(file="temp_figures/03/QE_series/03_venn_QE_S.pdf")
jco_col = pal_jco("default", alpha = 0.5)(10)[c(1,2,3,4,5)]
temp <- venn.diagram(
  x = list(
    QE = as.vector(qe.px.fot$GeneSymbol),
    QE_Plus = as.vector(qeplus.px.fot$GeneSymbol),
    QE_HF = as.vector(qehf.px.fot$GeneSymbol),
    Exploris480 = as.vector(e480.px.fot$GeneSymbol),
    QE_HFX = as.vector(qehfx.px.fot$GeneSymbol)
  ),
  # filename = 'temp_figures/03/QE_series/03_venn_QE_S.tiff',
  filename = NULL,
  main = 'QE-Series (G1->G5)',
  # main.just = c(0, 5),
  main.pos = c(0.5,1.05),
  main.cex = 0.4,
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
grid.draw(temp)
dev.off()

#2 identification frequency
QE.row_id = apply(qe.px.fot[,-1], 1, function(x){length(which(x>0))})
QE_Plus.row_id = apply(qeplus.px.fot[,-1], 1, function(x){length(which(x>0))})
QE_HF.row_id = apply(qehf.px.fot[,-1], 1, function(x){length(which(x>0))})
QE_HFX.row_id = apply(qehfx.px.fot[,-1], 1, function(x){length(which(x>0))})
Exploris480.row_id = apply(e480.px.fot[,-1], 1, function(x){length(which(x>0))})










#3 intensity distribution
get_expr_label <- function(x){
  # x = qe.row_mean
  x.quantile = as.vector(quantile(x, probs = c(1/3, 2/3)))
  x.l_index = which(x<x.quantile[1])
  x.h_index = which(x>x.quantile[2])
  x.expr_label = rep('Medium', length(x))
  x.expr_label[x.l_index] = 'Low'
  x.expr_label[x.h_index] = 'High'
  x.expr_label = factor(x.expr_label, levels = c('Low', 'Medium', 'High'))
  return(x.expr_label)
}


QE.row_mean = rowMeans(qe.px.fot[,-1])
QE.row_expr_label = get_expr_label(QE.row_mean)

QE_Plus.row_mean = rowMeans(qeplus.px.fot[,-1])
QE_Plus.row_expr_label = get_expr_label(QE_Plus.row_mean)

QE_HF.row_mean = rowMeans(qehf.px.fot[,-1])
QE_HF.row_expr_label = get_expr_label(QE_HF.row_mean)

QE_HFX.row_mean = rowMeans(qehfx.px.fot[,-1])
QE_HFX.row_expr_label = get_expr_label(QE_HFX.row_mean)

Exploris480.row_mean = rowMeans(e480.px.fot[,-1])
Exploris480.row_expr_label = get_expr_label(Exploris480.row_mean)

QE.df = data.frame(
  GeneSymbol = qe.px.fot$GeneSymbol,
  RowID = QE.row_id,
  ExprLab = QE.row_expr_label
)

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

QE_HFX.df = data.frame(
  GeneSymbol = qehfx.px.fot$GeneSymbol,
  RowID = QE_HFX.row_id,
  ExprLab = QE_HFX.row_expr_label
)

Exploris480.df = data.frame(
  GeneSymbol = e480.px.fot$GeneSymbol,
  RowID = Exploris480.row_id,
  ExprLab = Exploris480.row_expr_label
)


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
  
  
  # qe.px.col_gps, qe.df
  qe.repro = get_px_reproducibility(qe.px.col_gps, ref_df = NULL, ref_lab = NULL)
  qe.repro.h = get_px_reproducibility(qe.px.col_gps, ref_df = QE.df, ref_lab = 'High')
  qe.repro.m = get_px_reproducibility(qe.px.col_gps, ref_df = QE.df, ref_lab = 'Medium')
  qe.repro.l = get_px_reproducibility(qe.px.col_gps, ref_df = QE.df, ref_lab = 'Low')
  df.qe = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(qe.repro.l, qe.repro.m, qe.repro.h, qe.repro)
  )
  df.qe$Intensity = factor(df.qe$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  
  # qeplus.px.col_gps, qeplus.df
  qeplus.repro = get_px_reproducibility(qeplus.px.col_gps, ref_df = NULL, ref_lab = NULL)
  qeplus.repro.h = get_px_reproducibility(qeplus.px.col_gps, ref_df = QE_Plus.df, ref_lab = 'High')
  qeplus.repro.m = get_px_reproducibility(qeplus.px.col_gps, ref_df = QE_Plus.df, ref_lab = 'Medium')
  qeplus.repro.l = get_px_reproducibility(qeplus.px.col_gps, ref_df = QE_Plus.df, ref_lab = 'Low')
  df.qeplus = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(qeplus.repro.l, qeplus.repro.m, qeplus.repro.h, qeplus.repro)
  )
  df.qeplus$Intensity = factor(df.qeplus$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  
  
  # qehf.px.col_gps, qehf.df
  qehf.repro = get_px_reproducibility(qehf.px.col_gps, ref_df = NULL, ref_lab = NULL)
  qehf.repro.h = get_px_reproducibility(qehf.px.col_gps, ref_df = QE_HF.df, ref_lab = 'High')
  qehf.repro.m = get_px_reproducibility(qehf.px.col_gps, ref_df = QE_HF.df, ref_lab = 'Medium')
  qehf.repro.l = get_px_reproducibility(qehf.px.col_gps, ref_df = QE_HF.df, ref_lab = 'Low')
  df.qehf = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(qehf.repro.l, qehf.repro.m, qehf.repro.h, qehf.repro)
  )
  df.qehf$Intensity = factor(df.qehf$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  
  # qehfx.px.col_gps, qehfx.df
  qehfx.repro = get_px_reproducibility(qehfx.px.col_gps, ref_df = NULL, ref_lab = NULL)
  qehfx.repro.h = get_px_reproducibility(qehfx.px.col_gps, ref_df = QE_HFX.df, ref_lab = 'High')
  qehfx.repro.m = get_px_reproducibility(qehfx.px.col_gps, ref_df = QE_HFX.df, ref_lab = 'Medium')
  qehfx.repro.l = get_px_reproducibility(qehfx.px.col_gps, ref_df = QE_HFX.df, ref_lab = 'Low')
  df.qehfx = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(qehfx.repro.l, qehfx.repro.m, qehfx.repro.h, qehfx.repro)
  )
  df.qehfx$Intensity = factor(df.qehfx$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  
  
  # e480.px.col_gps, e480.df
  e480.repro = get_px_reproducibility(e480.px.col_gps, ref_df = NULL, ref_lab = NULL)
  e480.repro.h = get_px_reproducibility(e480.px.col_gps, ref_df = Exploris480.df, ref_lab = 'High')
  e480.repro.m = get_px_reproducibility(e480.px.col_gps, ref_df = Exploris480.df, ref_lab = 'Medium')
  e480.repro.l = get_px_reproducibility(e480.px.col_gps, ref_df = Exploris480.df, ref_lab = 'Low')
  df.e480 = data.frame(
    Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
    Reproducibity = c(e480.repro.l, e480.repro.m, e480.repro.h, e480.repro)
  )
  df.e480$Intensity = factor(df.e480$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))
  
  
  
  
  
  
  df = rbind.data.frame(
    df.qe, 
    df.qeplus,
    df.qehf,
    df.qehfx,
    df.e480
  )
  df$Instrument = rep(
    c('G1', 'G2', 'G3', 'G4', 'G5'), 
    each = 24
  )
  df$Instrument = factor(df$Instrument, levels = c('G1', 'G2', 'G3', 'G4', 'G5'))
  
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
  if(T){
    qe.px.fot.row_id = apply(qe.px.fot[,-1], 1, function(x){length(which(x>0))})
    qeplus.px.fot.row_id = apply(qeplus.px.fot[,-1], 1, function(x){length(which(x>0))})
    qehf.px.fot.row_id = apply(qehf.px.fot[,-1], 1, function(x){length(which(x>0))})
    qehfx.px.fot.row_id = apply(qehfx.px.fot[,-1], 1, function(x){length(which(x>0))})
    e480.px.fot.row_id = apply(e480.px.fot[,-1], 1, function(x){length(which(x>0))})
    px.fot = merge(
      qe.px.fot[qe.px.fot.row_id>=1,], 
      merge(
        qeplus.px.fot[qeplus.px.fot.row_id>=1,], 
        merge(
          qehf.px.fot[qehf.px.fot.row_id>=1,], 
          merge(
            qehfx.px.fot[qehfx.px.fot.row_id>=1,], 
            e480.px.fot[e480.px.fot.row_id>=1,],
            by = 'GeneSymbol', all = T
          ),
          by = 'GeneSymbol', all = T
        ),
        by = 'GeneSymbol', all = T
      ), 
      by = 'GeneSymbol', all = T
    )
    
    
    
    
    px.fot.mat = as.matrix(px.fot[,-1])
    colnames(px.fot.mat) = c(
      paste('G1', c('R1', 'R2', 'R3'), sep = '_'),
      paste('G2', c('R1', 'R2', 'R3'), sep = '_'),
      paste('G3', c('R1', 'R2', 'R3'), sep = '_'),
      paste('G4', c('R1', 'R2', 'R3'), sep = '_'),
      paste('G5', c('R1', 'R2', 'R3'), sep = '_')
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
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )
    p1 <- p1 + scale_fill_continuous(
      limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1),
      guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)
    ) + labs(fill = "PCorr")
    p1
    
    # col<- colorRampPalette(c("blue", "white", "red"))(20)
    # heatmap(x = cormat, col = col, symm = TRUE)
  }
}



library(ggforce)
# QE
df = QE.df
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
  title = paste('G1 (Total=', nrow(df), ')', sep = '')
)

# QE-Plus
df = QE_Plus.df
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
  title = paste('G2 (Total=', nrow(df), ')', sep = '')
)


# QE-HF
df = QE_HF.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
PieDonut_zdd(
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(3,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = paste('G3 (Total=', nrow(df), ')', sep = '')
)


# QE-HFX
df = QE_HFX.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
PieDonut_zdd(
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(3,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = paste('G4 (Total=', nrow(df), ')', sep = '')
)

# QE-HFX
df = Exploris480.df
df$IDF = paste('IDF=', df$RowID, sep = '')
df$D5 = df$ExprLab
PieDonut_zdd(
  df, 
  aes(pies=D5, donuts=IDF),
  ratioByGroup = TRUE,
  labelposition = 0,
  # explode = 1, 
  selected = c(3,6,9),
  explodeDonut = T,
  r0 = 0.2, r1 = 0.8, r2 = 1.2,
  title = paste('G5 (Total=', nrow(df), ')', sep = '')
)






generate_sankey_df <- function(df1, df2, instruments){
  # instruments = c('QE', 'QE-Plus')
  # df1 = QE.df
  # df2 = QE_Plus.df
  df.merge = merge(df1, df2, by = 'GeneSymbol', all = T)
  df.merge$ExprLab.x = as.vector(df.merge$ExprLab.x)
  df.merge$ExprLab.y = as.vector(df.merge$ExprLab.y)
  df.merge$RowID.x[which(is.na(df.merge$RowID.x))] = 'NA'
  df.merge$ExprLab.x[which(is.na(df.merge$ExprLab.x))] = 'NA'
  df.merge$RowID.y[which(is.na(df.merge$RowID.y))] = 'NA'
  df.merge$ExprLab.y[which(is.na(df.merge$ExprLab.y))] = 'NA'
  
  lable_list = c('High', 'Medium', 'Low', 'NA')
  labels_df = NULL
  for (i in 1:4) {
    df1.label = lable_list[i]
    for (j in 1:4) {
      df2.label = lable_list[j]
      i_j.index = which(
        df.merge$ExprLab.x == df1.label & df.merge$ExprLab.y == df2.label
      )
      i_j_labels = c(instruments[1], df1.label, instruments[2], df2.label, length(i_j.index))
      labels_df = rbind(labels_df,i_j_labels )
    }
  }
  labels_df = as.data.frame(labels_df[-16,])
  colnames(labels_df) = c('MS1', 'Lab1', 'MS2', 'Lab2', 'Value')
  return(labels_df)
}

construct_link_df <- function(df){
  # df = qe_qeplus.df
  df$source = paste(df$MS1, df$Lab1, sep = '-')
  df$target = paste(df$MS2, df$Lab2, sep = '-')
  df$group = df$source
  df$value = df$Value
  return(df)
}

qe_qeplus.df = construct_link_df(generate_sankey_df(QE.df, QE_Plus.df, c('G1', 'G2')))
qe_qehf.df = construct_link_df(generate_sankey_df(QE.df, QE_HF.df, c('G1', 'G3')))
qe_qehfx.df = construct_link_df(generate_sankey_df(QE.df, QE_HFX.df, c('G1', 'G4')))
qe_e480.df = construct_link_df(generate_sankey_df(QE.df, Exploris480.df, c('G1', 'G5')))

qeplus_qehf.df = construct_link_df(generate_sankey_df(QE_Plus.df, QE_HF.df, c('G2', 'G3')))
qeplus_qehfx.df = construct_link_df(generate_sankey_df(QE_Plus.df, QE_HFX.df, c('G2', 'G4')))
qeplus_e480.df = construct_link_df(generate_sankey_df(QE_Plus.df, Exploris480.df, c('G2', 'G5')))

qehf_qehfx.df = construct_link_df(generate_sankey_df(QE_HF.df, QE_HFX.df, c('G3', 'G4')))
qehf_e480.df = construct_link_df(generate_sankey_df(QE_HF.df, Exploris480.df, c('G3', 'G5')))

qehfx_e480.df = construct_link_df(generate_sankey_df(QE_HFX.df, Exploris480.df, c('G4', 'G5')))



if(T){
  
  links = qe_qeplus.df
  links = rbind.data.frame(
    qe_qeplus.df, 
    qeplus_qehf.df,
    qehf_qehfx.df,
    qehfx_e480.df
  )
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
  my_color <- 'd3.scaleOrdinal() .domain([
  "G1-High", "G1-Medium", "G1-Low", "G1-NA", 
  "G2-High", "G2-Medium", "G2-Low", "G2-NA",
  "G3-High", "G3-Medium", "G3-Low", "G3-NA",
  "G4-High", "G4-Medium", "G4-Low", "G4-NA",
  "G5-High", "G5-Medium", "G5-Low", "G5-NA",
  ]) .range([
  "#F9918A", "#81B0FF", "#33C860", "grey", 
  "#F9918A", "#81B0FF", "#33C860", "grey",
  "#F9918A", "#81B0FF", "#33C860", "grey",
  "#F9918A", "#81B0FF", "#33C860", "grey",
  "#F9918A", "#81B0FF", "#33C860", "grey"
  ])'
  p1 <- sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "IDsource", Target = "IDtarget",
    Value = "value", 
    # width = 800,
    NodeID = "name", 
    fontSize = 15,
    colourScale = my_color,
    LinkGroup="group",
    sinksRight = FALSE
  )
  p1
}

# chordDiagram

#### PCA #####
if(T){
  # PCA
  qe.fot = get_fot_df_by_expnames(fot_data, meta_data.NPS.QE$Exp_code)
  qeplus.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEPlus$Exp_code)
  qehf.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEHF$Exp_code)
  qehfx.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.QEHFX$Exp_code)
  e480.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.E480$Exp_code)
  
  qe.fot.ref = get_fot_ref_df(qe.fot)
  qeplus.fot.ref = get_fot_ref_df(qeplus.fot)
  qehf.fot.ref = get_fot_ref_df(qehf.fot)
  qehfx.fot.ref = get_fot_ref_df(qehfx.fot)
  e480.fot.ref = get_fot_ref_df(e480.fot)
  
  pdf(
    normalizePath(
      file.path(BASE_DIR, 'temp_figures', '03', 'QE_series', 'PCA_QE_S.pdf'), mustWork = F
    ),
    height = 8, width = 8
  )
  g1.snr_list = plot_single_pca(qe.fot, qe.fot.ref, main_lab = 'G1', p_flags, p_colors)
  g2.snr_list = plot_single_pca(qeplus.fot, qeplus.fot.ref, main_lab = 'G2', p_flags, p_colors)
  g3.snr_list = plot_single_pca(qehf.fot, qehf.fot.ref, main_lab = 'G3', p_flags, p_colors)
  g4.snr_list = plot_single_pca(qehfx.fot, qehfx.fot.ref, main_lab = 'G4', p_flags, p_colors)
  g5.snr_list = plot_single_pca(e480.fot, e480.fot.ref, main_lab = 'G5', p_flags, p_colors)
  dev.off()
  snr_list = rbind(
    g1.snr_list, g2.snr_list, g3.snr_list,
    g4.snr_list, g5.snr_list
  )
  
  
  
  
  get_fot_ref_df <- function(df){
    # qe.fot
    df.row_mean = rowMeans(df[,-1])
    df.row_expr_label = get_expr_label(df.row_mean)
    df.df = data.frame(
      GeneSymbol = df$GeneSymbol,
      RowID = apply(df[,-1], 1, function(x){length(which(x>0))}),
      MeanFOT = df.row_mean,
      ExprLab = df.row_expr_label
    )
    return(df.df)
  }
  
  # Single PCA
  # qe
  plot_single_pca <- function(fot_df, fot_df.ref, main_lab = 'G1', p_flags, p_colors, cex=2){
    # fot_df = qe.fot
    # fot_df.ref = qe.fot.df
    
    fot_df.mat = as.matrix(fot_df[,-1])
    fot_df.mat.scale = t(scale(t(fot_df.mat)))
    fot_group = rep(p_flags, each = 3)
    par(mfrow = c(2,2))
    snr_list = NULL
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
      main = paste(main_lab, '; Low-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main,
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
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
      main = paste(main_lab, '; Medium-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
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
      main = paste(main_lab, '; High-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
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
      main = paste(main_lab, '; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
    }
    names(snr_list) = c('Low', 'Medium', 'High', 'Global')
    return(snr_list)
  }
}  


#### CV ####
if(T){
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px = px_list[1]
  
  # QE
  qe.px.index = which(meta_data.NPS.QE$UPerson == px)
  qe.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPS.QE$Exp_code[qe.px.index])
  
  # QE-Plus
  qeplus.px.index = which(meta_data.NPB.QEPlus$UPerson == px)
  qeplus.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEPlus$Exp_code[qeplus.px.index])
  
  # QE-HF
  qehf.px.index = which(meta_data.NPB.QEHF$UPerson == px)
  qehf.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.QEHF$Exp_code[qehf.px.index])
  
  # QE-HFX
  qehfx.px.index = which(meta_data.TMO.QEHFX$UPerson == px)
  qehfx.px.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.QEHFX$Exp_code[qehfx.px.index])
  
  # Exploris480
  e480.px.index = which(meta_data.TMO.E480$UPerson == px)
  e480.px.fot = get_fot_df_by_expnames(fot_data, meta_data.TMO.E480$Exp_code[e480.px.index])
  
  get_cv_list_set <- function(df, min_id = 3){
    
    df.row_mean = rowMeans(df[,-1])
    df.row_id = apply(df[,-1], 1, function(x){length(which(x>0))})
    df.row_expr_label = get_expr_label(df.row_mean)
    df.ref_df = data.frame(
      GeneSymbol = df$GeneSymbol,
      RowID = df.row_id,
      ExprLab = df.row_expr_label
    )
    
    kept.index = which(df.ref_df$RowID>=min_id)
    df1 = df[kept.index,]
    df1.ref_df = df.ref_df[kept.index,]
    
    # CV
    h_index = which(df1.ref_df$ExprLab=='High')
    m_index = which(df1.ref_df$ExprLab=='Medium')
    l_index = which(df1.ref_df$ExprLab=='Low')
    
    
    get_cv_list <- function(df){
      # df = fusion.px.fot.common
      df.row_cv = apply(df[,-1], 1, function(x){
        sd(x)/mean(x)
      })
      df.row_cv = round(df.row_cv, 2)
      return(df.row_cv)
    }
    
    df1.g.cv = get_cv_list(df1)
    df1.h.cv = get_cv_list(df1[h_index,])
    df1.m.cv = get_cv_list(df1[m_index,])
    df1.l.cv = get_cv_list(df1[l_index,])
    
    cv_list = list(
      g_cv = df1.g.cv, h_cv = df1.h.cv,
      m_cv = df1.m.cv, l_cv = df1.l.cv
    )
    
    # lapply(cv_list, quantile)
    return(cv_list)
  }
  
  qe.px.fot.cv_list = get_cv_list_set(qe.px.fot, min_id = 3)
  qeplus.px.fot.cv_list = get_cv_list_set(qeplus.px.fot, min_id = 3)
  qehf.px.fot.cv_list = get_cv_list_set(qehf.px.fot, min_id = 3)
  qehfx.px.fot.cv_list = get_cv_list_set(qehf.px.fot, min_id = 3)
  e480.px.fot.cv_list = get_cv_list_set(e480.px.fot, min_id = 3)
  
  
  
  # '#33C860', # 绿色
  # '#81B0FF', # 蓝色
  # '#F9918A'  # 红色
  
  if(T){
    bp <- boxplot(
      NA, 
      xlim = c(0, 21.5), ylim = c(0, 2), ylab = 'Coefficient of Variation',
      xaxt = 'n', yaxt = 'n'
    )
    boxplot(
      qe.px.fot.cv_list$l_cv, qe.px.fot.cv_list$m_cv,
      qe.px.fot.cv_list$h_cv, qe.px.fot.cv_list$g_cv,
      at = c(0.25, 1.25, 2.25, 3.25), boxwex = 0.5,
      add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
    )
    
    boxplot(
      qeplus.px.fot.cv_list$l_cv, qeplus.px.fot.cv_list$m_cv,
      qeplus.px.fot.cv_list$h_cv, qeplus.px.fot.cv_list$g_cv,
      at = c(0.25, 1.25, 2.25, 3.25)+4.5, boxwex = 0.5,
      add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
    )
    
    boxplot(
      qehf.px.fot.cv_list$l_cv, qehf.px.fot.cv_list$m_cv,
      qehf.px.fot.cv_list$h_cv, qehf.px.fot.cv_list$g_cv,
      at = c(0.25, 1.25, 2.25, 3.25)+4.5*2, boxwex = 0.5,
      add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
    )
    
    boxplot(
      qehfx.px.fot.cv_list$l_cv, qehfx.px.fot.cv_list$m_cv,
      qehfx.px.fot.cv_list$h_cv, qehfx.px.fot.cv_list$g_cv,
      at = c(0.25, 1.25, 2.25, 3.25)+4.5*3, boxwex = 0.5,
      add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
    )
    
    boxplot(
      e480.px.fot.cv_list$l_cv, e480.px.fot.cv_list$m_cv,
      e480.px.fot.cv_list$h_cv, e480.px.fot.cv_list$g_cv,
      at = c(0.25, 1.25, 2.25, 3.25)+4.5*4, boxwex = 0.5,
      add = T, xaxt = 'n', yaxt = 'n', col = c('#33C860', '#81B0FF', '#F9918A', '#FFCC00')
    )
    
    
    
    
    
    abline(v = c(4, 8.5, 13, 17.5), lty = 3, lwd = 2, col = 'grey')
    axis(
      1, at = c(1.75, 6.25, 10.75, 15.25, 19.75), 
      labels = c('G1', 'G2', 'G3', 'G4', 'G5')
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
  # get_probs_table <- function(df){
  #   # qe_qeplus.df
  #   # qeplus_qehf.df
  #   # qehf_qehfx.df
  #   # qehfx_e480.df
  #   df = df[,c('source', 'target', 'value')]
  #   del.index = apply(df, 1, function(x){
  #     x = as.vector(unlist(x))
  #     if(grepl('NA', x[1])){
  #       return(T)
  #     }else{
  #       return(F)
  #     }
  #   })
  #   df1 = df[-which(del.index),]
  #   df1.mat = matrix(as.numeric(as.vector(df1$value)), nrow = 3, ncol = 4, byrow = T)
  #   rownames(df1.mat) = unique(df1$source)
  #   colnames(df1.mat) = unique(df1$target)
  #   df1.mat.rowsum = apply(df1.mat, 1, sum)
  #   df1.mat.probs_table = df1.mat/df1.mat.rowsum
  #   df1.mat.probs_table = round(df1.mat.probs_table, 2)
  #   return(df1.mat.probs_table)
  # }
  
  get_probs_table <- function(df){
    # qe_qeplus.df
    # qeplus_qehf.df
    # qehf_qehfx.df
    # qehfx_e480.df
    df = qe_qeplus.df
    df = df[,c('source', 'target', 'value')]

    del.index = apply(df, 1, function(x){
      x = as.vector(unlist(x))
      if(grepl('NA', x[1])){
        return(T)
      }else{
        return(F)
      }
    })
    df1 = df#[-which(del.index),]
    df1.mat = matrix(as.numeric(as.vector(df1$value)), nrow = 3, ncol = 4, byrow = T)
    rownames(df1.mat) = unique(df1$source)
    colnames(df1.mat) = unique(df1$target)
    df1.mat.rowsum = apply(df1.mat, 1, sum)
    df1.mat.probs_table = df1.mat/df1.mat.rowsum
    df1.mat.probs_table = round(df1.mat.probs_table, 2)
    return(df1.mat.probs_table)
  }
  
  qe_qeplus.df.table = get_probs_table(qe_qeplus.df)
  qeplus_qehf.df.table = get_probs_table(qeplus_qehf.df)
  qehf_qehfx.df.table = get_probs_table(qehf_qehfx.df)
  qehfx_e480.df.table = get_probs_table(qehfx_e480.df)
  
  output_path = normalizePath(path = file.path(
    BASE_DIR, 'sankey','QE-series-sankey-probs.xlsx'
  ),mustWork = F)
  
  write.xlsx(
    qe_qeplus.df.table, output_path, sheetName = 'qe_qeplus'
  )
  write.xlsx(
    qeplus_qehf.df.table, output_path, sheetName = 'qeplus_qehf', append = T
  )
  write.xlsx(
    qehf_qehfx.df.table, output_path, sheetName = 'qehf_qehfx', append = T
  )
  write.xlsx(
    qehfx_e480.df.table, output_path, sheetName = 'qehfx_e480', append = T
  )
  
}

