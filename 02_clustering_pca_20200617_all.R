library(factoextra)
library(dendextend)

meta_data.cp = meta_data[which(meta_data$Purpose == 'Replications'),]
meta_data.cp.uids = paste(
  meta_data.cp$Labcode, 
  meta_data.cp$SInstrument, 
  meta_data.cp$Batch, sep = '_'
)
uids = sort(unique(meta_data.cp.uids))
uids_count = length(uids)
uid.df.scale.list = list()
snr_list = NULL
for (i in 1:uids_count) {
  uid = uids[i] 
  uid.metadata = meta_data.cp[meta_data.cp.uids==uid,]
  # uid.metadata = uid.metadata[uid.metadata$Person!='P8',]
  # uid.metadata = uid.metadata[uid.metadata$Person=='P5' | uid.metadata$Person=='P6',]
  uid.metadata = uid.metadata[order(uid.metadata$Person),]
  uid.expnames = uid.metadata$Exp_code
  uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
  uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
  uid.df.mat = as.vector(uid.df[,-1])
  uid.df.mat.scale = t(scale(t(uid.df.mat)))
  
  mat = uid.df.mat.scale
  fot_group = uid.metadata$Person
  snr = round(snr_function(mat, fot_group),2)
  colors = rep(p_colors, each=3)
  pca <- prcomp(((t(mat))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste(uid, '; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
  plot(
    pca_predict.2d, 
    t = 'n', 
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d, pch = 16, col = colors)
  print(i)
  snr_list = c(snr_list, snr)
}
names(snr_list) = uids
snr_list.asec.order = order(snr_list)

snr_list = sort(snr_list)

df = data.frame(
  UID = uids[snr_list.asec.order], 
  SNR = snr_list
)
# df = data.frame(
#   UID = c(uids, 'FDU_QE-HFX_3', 'FDU_QE-HFX_4'),
#   SNR = c(snr_list, 20.3, 16.3)
# )
df = df[order(df$SNR),]
df$SNR = as.numeric(df$SNR)

jco_col = pal_jco("default", alpha = 0.5)(10)#[1:5]
raw_cols_rgb = col2rgb('firebrick')
alpha_cols = rgb(
  raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], 
  raw_cols_rgb[3L, ], alpha = 0.7 * 255L, names = 'zdd_firebrick', 
  maxColorValue = 255L
)
jco_col.list = c(
  rep(jco_col[3], 5),
  rep(jco_col[5], 6),
  rep(jco_col[1], 9),
  rep(jco_col[4], 3),
  rep(jco_col[9], 1)
)

jco_col.list = c(
  rep(jco_col[3], 5),
  rep(jco_col[5], 6),
  rep(jco_col[1], 9),
  rep(jco_col[4], 4),
  rep(alpha_cols, 2)
)

jco_col.list = c(
  rep(jco_col[3], 5),
  rep(jco_col[5], 6),
  rep(jco_col[1], 9),
  rep(jco_col[4], 3),
  alpha_cols
)

df$Level =  c(
  rep("D", 5),
  rep('C', 6),
  rep("B", 9),
  rep('A', 3),
  rep('A+', 1)
)
df$Level = factor(
  df$Level, c('D', 'C', 'B', 'A', 'A+')
)

Company =  apply(df, 1, function(x){
  strsplit(x, '_')[[1]][1]
})
Instrument = apply(df, 1, function(x){
  strsplit(x, '_')[[1]][2]
})
Batch = apply(df, 1, function(x){
  strsplit(x, '_')[[1]][3]
})

df1 = data.frame(
  df, Company, Instrument, Batch
)

df$UID = paste('T', seq(1,24), sep = '')
ggdotchart(
  df, x = "UID", y = "SNR",
  color = jco_col.list,                                # Color by groups
  palette = jco_col.list, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
  # sorting = "ascending",                        # Sort value in descending order
  add = "segments",                             # Add segments from y = 0 to dots
  ggtheme = theme_pubr(),                        # ggplot2 theme
  # group = "Level",                                # Order by groups
  dot.size = 6,                                 # Large dot size
  label = round(df$SNR,1),
  add.params = list(color = jco_col.list, size = 1.5), # Change segment color and size,
  xlab = '',
  ylab = 'Signal-to-noise Ratio',
  ylim = c(0, 30),
  font.label = list(color = "black", size = 10, vjust = -1)
)

write.csv(
  df1, '26_SNR_Figures.csv', row.names = F
)



pdf('24_all_pca.pdf', height = 5, width = 5)

for (i in 1:uids_count) {
  uid = uids[snr_list.asec.order[i]] 
  uid.metadata = meta_data.cp[meta_data.cp.uids==uid,]
  uid.metadata = uid.metadata[order(uid.metadata$Person),]
  uid.expnames = uid.metadata$Exp_code
  uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
  uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
  uid.df.mat = as.vector(uid.df[,-1])
  uid.df.mat.scale = t(scale(t(uid.df.mat)))
  
  mat = uid.df.mat.scale
  fot_group = uid.metadata$Person
  snr = round(snr_function(mat, fot_group),2)
  colors = rep(p_colors, each=3)
  pca <- prcomp(((t(mat))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste(uid, '; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
  plot(
    pca_predict.2d, 
    t = 'n', 
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d, pch = 16, col = colors)
  print(i)
}
dev.off()






#### FDU new ####
if(T){
  fdu_new.fot_path = normalizePath(
    file.path(INPUT_DIR,'maqc_proteome_FDU_new_FOT_us=1.xlsx'), mustWork = F
  )
  metadata_path =  normalizePath(
    file.path(INPUT_DIR,'maqc_proteome_metadata_uniform.xlsx'), mustWork = F
  )
  
  fdu_new.fot_data = read_excel(fdu_new.fot_path)
  fdu_new.meta_data = read_excel(metadata_path, sheet = 'FDU_HFX_20200617')
  
  if(T){

    uid.metadata = fdu_new.meta_data
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fdu_new.fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    
    mat = uid.df.mat.scale
    fot_group = uid.metadata$Person
    fdu_new.snr = round(snr_function(mat, fot_group),2)
    colors = rep(p_colors, each=3)
    pca <- prcomp(((t(mat))), center = F, scale = F)
    importance <- summary(pca)$importance
    pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
    pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
    pca_predict <- predict(pca)
    pca_predict.2d = pca_predict[,c(1,2)]
    main = paste('FDU_QE-HFX_2', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
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
  print(fdu_new.snr)
  names(fdu_new.snr) = 'FDU_QE-HFX_2'
  
}




df1 = rbind(
  df, c('FDU_QE-HFX_2', fdu_new.snr, 'A')
)

jco_col = pal_jco("default", alpha = 0.5)(10)#[1:5]
jco_col.list = c(
  rep(jco_col[3], 5),
  rep(jco_col[5], 6),
  rep(jco_col[1], 9),
  rep(jco_col[4], 3),
  rep(jco_col[9], 1)
)
ggdotchart(
  df1, x = "UID", y = "SNR",
  color = jco_col.list,                                # Color by groups
  palette = jco_col.list, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
  # sorting = "ascending",                        # Sort value in descending order
  add = "segments",                             # Add segments from y = 0 to dots
  ggtheme = theme_pubr(),                        # ggplot2 theme
  # group = "Level",                                # Order by groups
  dot.size = 6,                                 # Large dot size
  label = round(df$SNR),
  add.params = list(color = jco_col.list, size = 2), # Change segment color and size,
  xlab = '',
  ylab = 'Protein identification',
  ylim = c(0, 30),
  font.label = list(color = "black", size = 10, vjust = -1)
)







#### FDU new 2 ####
if(T){
  fdu_new.fot_path = normalizePath(
    file.path(INPUT_DIR,'maqc_proteome_FDU_new_4_FOT_us=1.xlsx'), mustWork = F
  )
  metadata_path =  normalizePath(
    file.path(INPUT_DIR,'maqc_proteome_metadata_uniform.xlsx'), mustWork = F
  )
  
  fdu_new.fot_data = read_excel(fdu_new.fot_path)
  fdu_new.meta_data = read_excel(metadata_path, sheet = 'FDU_HFX_20200710')
  
  write.csv(
    fdu_new.fot_data, 'maqc_proteome_data_for_ly1.csv', row.names = F
  )
  write.csv(
    fdu_new.meta_data, 'maqc_proteome_data_for_ly2.csv', row.names = F
  )
  
  
  if(T){
    
    uid.metadata = fdu_new.meta_data[fdu_new.meta_data$Order!=1,]
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fdu_new.fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.col_id = apply(uid.df.mat, 2, function(x){length(which(x>0))})
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    a = data.frame(
      uid.expnames, uid.df.mat.col_id
    )
    write.csv(
      a, 'a.csv', row.names = F
    )
    
    mat = uid.df.mat.scale
    fot_group = uid.metadata$Person
    fdu_new.snr = round(snr_function(mat, fot_group),2)
    colors = rep(p_colors, each=3)
    pca <- prcomp(((t(mat))), center = F, scale = F)
    importance <- summary(pca)$importance
    pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
    pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
    pca_predict <- predict(pca)
    pca_predict.2d = pca_predict[,c(1,2)]
    main = paste('FDU_QE-HFX_2', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    plot(
      pca_predict.2d, 
      t = 'n', 
      main = main, 
      xlab = pc1, ylab = pc2, 
      xlim = get_xlim(pca_predict.2d), 
      ylim = get_ylim(pca_predict.2d),
    )
    points(pca_predict.2d, pch = 16, col = colors)
    points(pca_predict.2d, pch = uid.metadata$Order+10, col = colors, cex =2)
    
    legend(
      'topleft', 
      legend = sort(unique(uid.metadata$Order)), 
      title = 'Order',
      pch = sort(unique(uid.metadata$Order))+10
    )
  }
  print(fdu_new.snr)
  names(fdu_new.snr) = 'FDU_QE-HFX_3'
  
}

text(pca_predict.2d, labels = rownames(pca_predict.2d), cex = 0.5, pos = 3)

df = data.frame(
  x = pca_predict.2d[,1],
  y = pca_predict.2d[,2],
  label = rownames(pca_predict.2d),
  Sample = factor(uid.metadata$Person),
  Order = uid.metadata$Order[match(rownames(df), uid.metadata$Exp_code)]-1
)

library(ggrepel)
main = paste('FDU_QE-HFX_4', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
p <- ggplot(df, aes(x=x, y=y))
p <- p + geom_point(aes(color = Sample, shape = factor(Order)), size=3, position = "identity")
p <- p + scale_color_manual(values=p_colors)
# p <- p + geom_text_repel(
#   aes(x, y, label= label), 
#   # fontface="bold", 
#   color="black", 
#   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
#   segment.colour = "grey50", size = 4
# ) 
p <- p + theme_classic(base_size = 16)
p <- p + xlim(-80, 60) + ylim(-60, 60) + xlab(pc1) + ylab(pc2) 
# p <- p + ggtitle(paste('SNR=',fdu_new.snr, sep = ''))
p + ggtitle(main) +  theme(plot.title = element_text(size = 15, face = "bold"))




# scale_color_discrete(name="cyl")
# scale_x_continuous(expand = c(0.5, 0))
# scale_y_continuous(expand = c(0.25, 0))



df = data.frame(
  UID = c(uids, 'FDU_QE-HFX_3'),
  SNR = c(snr_list, fdu_new.snr)
)
df = df[order(df$SNR),]
df$SNR = as.numeric(df$SNR)

jco_col = pal_jco("default", alpha = 0.6)(10)#[1:5]
jco_col.list = c(
  rep(jco_col[3], 4),
  rep(jco_col[5], 7),
  rep(jco_col[1], 9),
  rep(jco_col[4], 3),
  rep(jco_col[9], 2)
)

ggdotchart(
  df, x = "UID", y = "SNR",
  color = jco_col.list,                                # Color by groups
  palette = jco_col.list, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
  # sorting = "ascending",                        # Sort value in descending order
  add = "segments",                             # Add segments from y = 0 to dots
  ggtheme = theme_pubr(),                        # ggplot2 theme
  # group = "Level",                                # Order by groups
  dot.size = 5,                                 # Large dot size
  label = round(df$SNR,1),
  add.params = list(color = jco_col.list, size = 2), # Change segment color and size,
  xlab = '',
  ylab = 'SNR',
  ylim = c(0, 30),
  font.label = list(color = "black", size = 10, vjust = -1)
)


if(T){
  df = read_excel('24_SNR_Figures.xlsx')
  df = df[!is.na(df$Code),]
  df = df[order(df$Code),]
  df$SNR = as.numeric(df$SNR)
  
  jco_col = pal_jco("default", alpha = 0.6)(10)#[1:5]
  jco_col.list = c(
    rep(jco_col[1], 2),
    rep(jco_col[2], 5),
    rep(jco_col[4], 2)
  )
  jco_col.list = rep(jco_col[1], nrow(df))
  
  ggdotchart(
    df, x = "UID", y = "SNR",
    color = jco_col.list,                                # Color by groups
    palette = jco_col.list, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
    sorting = 'none',                        # Sort value in descending order
    add = "segments",                             # Add segments from y = 0 to dots
    ggtheme = theme_pubr(),                        # ggplot2 theme
    # group = "Level",                                # Order by groups
    dot.size = 5,                                 # Large dot size
    label = round(df$SNR,1),
    add.params = list(color = jco_col.list, size = 2), # Change segment color and size,
    order = df$UID,
    xlab = '',
    ylab = 'SNR',
    ylim = c(0, 30),
    font.label = list(color = "black", size = 10, vjust = -1)
  )
  
  
}



if(T){
  library(data.table)
  exprMat = mat
  group = fot_group
  IDs<- colnames(exprMat)
  IDs.group.mat<-data.table(
    IDs=IDs,
    group=group) 
  
  pca_prcomp <- prcomp(t(exprMat),retx=T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs,each=length(IDs)),
                        ID.B = rep(IDs,time=length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]
  
  dt.dist[,Dist:=sqrt(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  
  return(signoise)
}